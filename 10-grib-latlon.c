#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>

#define BUF8(X)  (buf[(X-5)])
#define BUF16(X) ((buf[(X-5)] << 8) | (buf[(X-4)]))
#define BUF32(X) ((buf[(X-5)] << 24) | (buf[(X-4)] << 16) | (buf[(X-3)] << 8) | (buf[(X-2)]))

int fill_clut(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a);
int hexdump(unsigned char *buf, int len);

int main(int argc, char **argv)
{
  FILE *fi, *fo;
  unsigned char buf[256000];
  unsigned char rawdata[256000];
  int tlen,stat;
  int refy,refm,refd,refh,refn,refs;
  int slen,dlen,grdw,grdh,lat0,lon0,lat1,lon1,dlat,dlon;
  int pnum,type, ft;
  int nbits,V,M;
  int ibuf,n;
  unsigned char clut_r[256],clut_g[256],clut_b[256],clut_a[256];


  argv++;
  while (--argc) {
    printf("%s\n",*argv);
    if ((fi = fopen(*argv,"r")) != NULL) {
      fread(buf,1,16,fi);	/* section 0 */
      tlen = BUF32(17);

      fread(&slen,4,1,fi);	/* section 1 */
      slen = ntohl(slen);
      fread(buf,1,slen-4,fi);

      refy = BUF16(13);
      refm = BUF8(15);
      refd = BUF8(16);
      refh = BUF8(17);
      refn = BUF8(18);
      refs = BUF8(19);

      stat = BUF8(20);
      if (stat == 1){
        break;
      }

      n = 0;
      while (fread(&slen,4,1,fi) > 0) {
        n++;
        if (slen == 0x37373737) {
          break;
        }
	slen = ntohl(slen);	/* section 3 */
        fread(buf,1,slen-4,fi);

        dlen = BUF32(7);
        grdw = BUF32(31);
        grdh = BUF32(35);
        // printf("grid number %d x %d = %d\n",grdw,grdh,dlen);
        lat0 = BUF32(47);
        lon0 = BUF32(51);
        lat1 = BUF32(56);
        lon1 = BUF32(60);
        dlon = BUF32(64);
        dlat = BUF32(68);


        fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 4 */
        fread(buf,1,slen-4,fi);
        pnum = BUF8(11);
        type = BUF8(12);
        ft = BUF32(19);
        if (ft < 0) {
          ft = 0x80000000 - ft;
        }
        //else {exit(0);}

        if (pnum == 203) {	// parameter number : 203 = prec level
          printf("%d, %d, %d, %d\n",lat0,lon0,lat1,lon1);
        }




        fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 5 */
        fread(buf,1,slen-4,fi);
        //printf("section 5 = %d, %d\n",slen-4,ft);
        nbits = BUF8(12);
        if (nbits != 8) {
          fprintf(stderr, "error nbits != 8\n");
          exit(1);
        }
        V = BUF16(13);
        M = BUF16(15);
        fill_clut(buf,clut_r,clut_g,clut_b,clut_a);

        fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 6 */
        fread(buf,1,slen-4,fi);
        //printf("V = %x (%d) ; M = %x (%d)\n", V,V,M,M);

        fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 7 */
        fread(buf,1,slen-4,fi);
        // hexdump(buf,slen-4);

      }

      fclose(fi);
    }
    argv++;
  }
}


int hexdump(unsigned char *buf, int len){
  int i;
  for(i = 0; i < len; i++) {
    if ((i % 16) == 0) {
      printf("%06X :", i);
    }
    else if ((i % 8) == 0) {
      printf(" ");
    }
    printf(" %02X",buf[i] & 0xff);
    if ((i % 16) == 15) {
      printf("\n");
    }
  }
  printf("\n");
  return 0;
}


int fill_clut(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a){
  int i, m, M, ibuf;
  M = BUF16(15);
  m = 0;
  clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
  for (i = 0; i < M; i++) {
    m = i + 1;
    ibuf = BUF16(18 + i*2);
    if (ibuf == 0) {	// 0 mm
      clut_r[m] = 192; clut_g[m] = 192; clut_b[m] = 192; clut_a[m] = 192;
    }
    else if (ibuf < 100) {	// 0 .. 1 mm
      clut_r[m] = 204; clut_g[m] = 255; clut_b[m] = 255; clut_a[m] = 255;
    }
    else if (ibuf < 400) {	// 1 .. 4 mm
      clut_r[m] = 0; clut_g[m] = 216; clut_b[m] = 235; clut_a[m] = 255;
    }
    else if (ibuf < 1600) {	// 4 .. 16 mm
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 235; clut_a[m] = 255;
    }
    else if (ibuf < 3200) {	// 16 .. 32 mm
      clut_r[m] = 235; clut_g[m] = 235; clut_b[m] = 0; clut_a[m] = 255;
    }
    else if (ibuf < 20000) {	// 32 .. mm
      clut_r[m] = 235; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 255;
    }
    else {
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
    }
  }
  return 0;
}
