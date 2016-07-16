/* grib-tile.c
 * input  : JMA high-resolution precipitation nowcast GRIB2 data 
 * output : single/large PGM or PNG image
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#undef USE_PNG	/* Now use PGM due to process speed */
#ifdef USE_PNG
#include "png.h"
#endif

#define BUF8(X)  (buf[(X-5)])
#define BUF16(X) ((buf[(X-5)] << 8) | (buf[(X-4)]))
#define BUF32(X) ((buf[(X-5)] << 24) | (buf[(X-4)] << 16) | (buf[(X-3)] << 8) | (buf[(X-2)]))




char mypath[1024];
char basedir[1024];
unsigned char clut_r[256],clut_g[256],clut_b[256],clut_a[256];
int verbose = 0;



int hexdump(unsigned char *buf, int len)
{
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

int fill_clut_mono(unsigned char *buf,unsigned char *clut_a)
{
  /* replace to simple conversion
   * 0,0.5,1.0,1.5,...    V = i * 2;  
   */

  int i, m, l, M, ibuf;
  M = BUF16(15);
  m = 0;
  clut_a[m] = 0;
  for (i = 0; i < M; i++) {
    m = i + 1;
    ibuf = BUF16(18 + i*2);
    l = ibuf * 2 / 100;
    if (l > 255) {
      l = 255;
    }
    clut_a[m] = l;
  }
  return 0;
}

int fill_clut(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a)
{
  int i, m, M, ibuf;
  M = BUF16(15);
  m = 0;
  clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
  for (i = 0; i < M; i++) {
    m = i + 1;
    ibuf = BUF16(18 + i*2);
    if (ibuf < 30) {	// 0 mm
      clut_r[m] = 192; clut_g[m] = 192; clut_b[m] = 192; clut_a[m] = 0;
    }
    else if (ibuf < 100) {	// 0 .. 1 mm
      clut_r[m] = 250; clut_g[m] = 250; clut_b[m] = 250; clut_a[m] = 158;
    }
    else if (ibuf < 400) {	// 1 .. 4 mm
      clut_r[m] = 0; clut_g[m] = 216; clut_b[m] = 235; clut_a[m] = 158;
    }
    else if (ibuf < 1600) {	// 4 .. 16 mm
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 235; clut_a[m] = 158;
    }
    else if (ibuf < 3200) {	// 16 .. 32 mm
      clut_r[m] = 235; clut_g[m] = 235; clut_b[m] = 0; clut_a[m] = 158;
    }
    else if (ibuf < 30000) {	// 32 .. mm
      clut_r[m] = 235; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 158;
    }
    else {
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
    }
  }
  //m = 1;
  //clut_a[m] = 0;
  return 0;
}

/*
 * extract run-length
 */
int extract(unsigned char *src_buf, int src_len, int maxv, int nbits, unsigned char *dst_buf)
{
  int lngu;
  int lastv;
  int pos = 0;
  int i;
  int lngv,lngn;

  lngu = (1 << nbits) - 1 - maxv;
  /*
  if (verbose) {
    printf("src_len = %d, maxv = %d, nbits = %d, lngu = %d\n", src_len, maxv, nbits, lngu);
  }
  */

  lngv = 0;
  for(i = 0; i < src_len; i++) {
    if (src_buf[i] <= maxv) {
      while (lngv-- > 0) {
        dst_buf[pos++] = lastv;
      }
      lastv = dst_buf[pos++] = src_buf[i];
      lngn = 1;
      lngv = 0;
    }
    else {
      lngv = lngv + (lngn * (src_buf[i] - (maxv + 1)));
      lngn *= lngu;
    }
  }
  while (lngv-- > 0) {
    dst_buf[pos++] = lastv;
  }
  return pos;
}

#ifdef USE_PNG
int png_save(char *filename, unsigned char *imgbuf, int total_gridw, int total_gridh)
{
  FILE *fo;
  png_structp     png_ptr;
  png_infop       info_ptr;
  png_bytep *rowp;
  int j;

  if ((fo = fopen(filename,"w")) == NULL) {
    perror(filename);
    return -1;
  }

  rowp = (png_bytepp)malloc(total_gridh * sizeof(png_bytep));
  for (j = 0; j < total_gridh; j++) {
    rowp[j] = imgbuf + (total_gridw * j) * 4;
  }

  if ((png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL)) == NULL){
    fprintf(stderr,"cannot allocate png_structp ");
    return -1;
  }
  if ((info_ptr = png_create_info_struct(png_ptr)) == NULL) {
    fprintf(stderr,"cannot allocate png_infop ");
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    return -1;
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr,  &info_ptr);
    fclose(fo);
    return -1;
  }
  png_init_io(png_ptr, fo);
  png_set_IHDR(png_ptr, info_ptr, total_gridw, total_gridh,
               8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_set_bgr(png_ptr);
  png_set_swap_alpha(png_ptr);
  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, rowp);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fo);
  free(rowp);
  return 0;
}
#endif

int pgm_save(char *filename, unsigned char *imgbuf, int total_gridw, int total_gridh)
{
  FILE *fo;

  if ((fo = fopen(filename,"w")) == NULL) {
    perror(filename);
    return -1;
  }
  fprintf(fo,"P5\n%d %d\n255\n",total_gridw, total_gridh);
  fwrite(imgbuf,1,total_gridw*total_gridh,fo);
  fclose(fo);
  return 0;
}

int grib_read(char *filename)
{
  FILE *fi, *fo, *fm;
  unsigned char buf[256000];
  unsigned char rawdata[256000];
  unsigned char *imgbuf, *p, *q, r,g,b,a;
  char msg[128];
  char fn[1024];
  int stat;
  struct tm reftm;
  time_t reftime,tval;
  int slen,dlen,grdw,grdh,lat0,lon0,lat1,lon1,dlat,dlon;
  int lat_min, lat_max, lon_min, lon_max, dlat_min, dlon_min;
  int total_gridw, total_gridh;
  int pnum,type,ft,current_ft;
  int nbits,V,M;
  long fpos;
  long lat,lon,i,j,px,py,qx,qy,ii,jj;

  if ((fi = fopen(filename,"r")) == NULL) {
    perror(filename);
    exit(1);
  }

  fread(buf,1,16,fi);	/* section 0 */
  if (strncmp(buf,"GRIB",4) != 0){
    fprintf(stderr,"%s is not GRIB file\n",filename);
    exit(1);
  }

  fread(&slen,4,1,fi);	/* section 1 */
  slen = ntohl(slen);
  fread(buf,1,slen-4,fi);
  reftm.tm_year = BUF16(13) - 1900;
  reftm.tm_mon  = BUF8(15) - 1;
  reftm.tm_mday = BUF8(16);
  reftm.tm_hour = BUF8(17);
  reftm.tm_min  = BUF8(18);
  reftm.tm_sec  = BUF8(19);
  reftime = timegm(&reftm);

  stat = BUF8(20);
  if (stat == 1){
    fprintf(stderr, "%s is test product. process canceled.\n", filename);
    exit(1);
  }

  fpos = ftell(fi);

  /* first, get lat,lon area */
  lat_max = -999999999;
  lon_max = -999999999;
  lat_min =  999999999;
  lon_min =  999999999;
  dlat_min = 999999999;
  dlon_min = 999999999;
  current_ft = -999;
  while (fread(&slen,4,1,fi) > 0) {
    if (slen == 0x37373737) {	// section 8 = termination
      break;
    }

    slen = ntohl(slen);	/* section 3 */
    fread(buf,1,slen-4,fi);
    dlen = BUF32(7);
    grdw = BUF32(31);
    grdh = BUF32(35);
    // if (verbose){printf("grid number %d x %d = %d\n",grdw,grdh,dlen);}
    lat0 = BUF32(47);
    lon0 = BUF32(51);
    lat1 = BUF32(56);
    lon1 = BUF32(60);
    dlon = BUF32(64);
    dlat = BUF32(68);
    // if (verbose){printf("%d,%d - %d,%d / %d,%d\n",lat0,lon0,lat1,lon1,dlat,dlon);}
    if (lat_max < lat0) {
      lat_max = lat0 + dlat / 2;
    }
    if (lat_min > lat1) {
      lat_min = lat1 - dlat / 2;
    }
    if (lon_min > lon0) {
      lon_min = lon0 - dlon / 2;
    }
    if (lon_max < lon1) {
      lon_max = lon1 + dlon / 2;
    }
    if (dlat_min > dlat) {
      dlat_min = dlat;
    }
    if (dlon_min > dlon) {
      dlon_min = dlon;
    }

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 4 */
    fread(buf,1,slen-4,fi);
    ft = BUF32(19);
    if (ft < 0) {
      ft = 0x80000000 - ft;
    }
    if (ft >= 0){	// only check observation
      break;
    }

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 5 */
    fread(buf,1,slen-4,fi);
    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 6 */
    fread(buf,1,slen-4,fi);
    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 7 */
    fread(buf,1,slen-4,fi);
  }

  total_gridw = (lon_max - lon_min) / dlon_min;
  total_gridh = (lat_max - lat_min) / dlat_min;

 
  if ((imgbuf = (unsigned char *)malloc(total_gridw * total_gridh * 4)) == NULL){
    perror("malloc");
    exit(1);
  }
  memset(imgbuf, 0x00, total_gridw * total_gridh);

  /* read again */
  fseek(fi,fpos,SEEK_SET);
  current_ft = -999;
  while (fread(&slen,4,1,fi) > 0) {
    if (slen == 0x37373737) {	// section 8 = termination
      tval = reftime + current_ft * 60;
      strftime(msg,sizeof(msg),"%Y-%m-%d-%H-%M",gmtime(&tval));
#ifdef USE_PNG
      sprintf(fn,"%s/%s.png",basedir,msg);
      if (png_save(fn,imgbuf,total_gridw,total_gridh) < 0) {
	return -1;
      }
#else
      sprintf(fn,"%s/%s.pgm",basedir,msg);
      if (pgm_save(fn,imgbuf,total_gridw,total_gridh) < 0) {
	return -1;
      }
#endif
      if (verbose){printf("%s\n",fn);}
      break;
    }

    slen = ntohl(slen);	/* section 3 */
    fread(buf,1,slen-4,fi);
    dlen = BUF32(7);
    grdw = BUF32(31);
    grdh = BUF32(35);
    lat0 = BUF32(47);
    lon0 = BUF32(51);
    lat1 = BUF32(56);
    lon1 = BUF32(60);
    dlon = BUF32(64);
    dlat = BUF32(68);

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 4 */
    fread(buf,1,slen-4,fi);
    pnum = BUF8(11);
    // parameter number : 203 = prec level
    type = BUF8(12);
    ft = BUF32(19);
    if (ft < 0) {
      ft = 0x80000000 - ft;
    }
    ft += 5;	// adjust for obervation end time
    if (ft != current_ft){
      if (current_ft >= 0) {
	tval = reftime + current_ft * 60;
	strftime(msg,sizeof(msg),"%Y-%m-%d-%H-%M",gmtime(&tval));
#ifdef USE_PNG
	sprintf(fn,"%s/%s.png",basedir,msg);
	if (png_save(fn,imgbuf,total_gridw,total_gridh) < 0) {
	  return -1;
	}
#else
	sprintf(fn,"%s/%s.pgm",basedir,msg);
	if (pgm_save(fn,imgbuf,total_gridw,total_gridh) < 0) {
	  return -1;
	}
#endif
	if (verbose){printf("%s\n",fn);}
      }
      memset(imgbuf, 0x00, total_gridw * total_gridh);
      current_ft = ft;
    }

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 5 */
    fread(buf,1,slen-4,fi);
    nbits = BUF8(12);
    if (nbits != 8) {
      fprintf(stderr, "error nbits %d != 8 not supported\n",nbits);
      exit(1);
    }
    V = BUF16(13);
    M = BUF16(15);
#ifdef USE_PNG
    fill_clut(buf,clut_r,clut_g,clut_b,clut_a);
#else
    fill_clut_mono(buf,clut_a);
#endif

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 6 */
    fread(buf,1,slen-4,fi);

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 7 */
    fread(buf,1,slen-4,fi);

    if (pnum == 214) {      // parameter number : 214 = error info
      continue;
    }
    if (extract(&buf[1], slen-5, V, nbits, rawdata) != dlen) {
      fprintf(stderr, "Error: unmatch data length %d\n", dlen);
      exit(1);
    }

    p = rawdata;
    lat = lat0 + dlat / 2;
    for (j = 0; j < grdh; j++) {
      py = total_gridh * (lat_max - lat) / (lat_max - lat_min);
      lat -= dlat;
      qy = total_gridh * (lat_max - lat) / (lat_max - lat_min);
      if (qy > total_gridh) {
        qy = total_gridh;
      }
      lon = lon0 - dlon / 2;
      for (i = 0; i < grdw; i++) {
        px = total_gridw * (lon - lon_min) / (lon_max - lon_min);
        lon += dlon;
        qx = total_gridw * (lon - lon_min) / (lon_max - lon_min);
        if (qx > total_gridw) {
          qx = total_gridw;
        }

#ifdef USE_PNG
        r = clut_r[*p];
        g = clut_g[*p];
        b = clut_b[*p];
        a = clut_a[*p];
        for (jj = py; jj < qy; jj++) {
          q = imgbuf + ((total_gridw * jj) + px) * 4;
          for (ii = px; ii < qx; ii++) {
            *q++ = a;
            *q++ = b;
            *q++ = g;
            *q++ = r;
          }
        }
#else
	a = clut_a[*p];
        for (jj = py; jj < qy; jj++) {
          q = imgbuf + ((total_gridw * jj) + px);
          for (ii = px; ii < qx; ii++) {
            *q++ = a;
          }
        }
#endif
        p++;
      }
    }
  }

  fclose(fi);

  sprintf(fn,"%s/info.txt",basedir);
  if ((fm = fopen(fn,"w")) == NULL) {
    perror(fn);
    exit(1);
  }

  fprintf(fm,"reftm\t%d\n",reftime);
  strftime(msg,sizeof(msg),"%Y-%m-%d-%H-%M",&reftm);
  fprintf(fm,"reftstr\t%s\n",msg);

  fprintf(fm,"wlon\t%d\n",lon_min);
  fprintf(fm,"elon\t%d\n",lon_max);
  fprintf(fm,"nlat\t%d\n",lat_max);
  fprintf(fm,"slat\t%d\n",lat_min);
  fprintf(fm,"gridw\t%d\n",total_gridw);
  fprintf(fm,"gridh\t%d\n",total_gridh);

  fclose(fm);
  return 0;
}

int main(int argc, char **argv)
{
  char *p;
  strcpy(mypath,*argv);
  if ((p = rindex(mypath,'/')) != NULL) {
    *p = '\0';
  }

  strcpy(basedir,".");
  int opt_basedir = 0;
  while (--argc) {
    argv++;
    if (opt_basedir == 1) {
      opt_basedir = 0;
      strcpy(basedir,*argv);
      continue;
    }
    else if (strncmp(*argv,"-d=",3) == 0) {
      strcpy(basedir,(*argv + 3));
      continue;
    }
    else if (strncmp(*argv,"-d",2) == 0) {
      opt_basedir = 1;
      continue;
    }
    else if (strncmp(*argv,"-v",2) == 0) {
      verbose = 1;
      continue;
    }
    else {
      grib_read(*argv);
    }
  }
}



