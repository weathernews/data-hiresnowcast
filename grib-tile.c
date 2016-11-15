/* grib-tile.c
 * input  : JMA high-resolution precipitation nowcast GRIB2 data 
 * output : GoogleMaps tile set
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <netinet/in.h>
#include "png.h"
#include "tile_info.h"

#define BUF8(X)  (buf[(X-5)])
#define BUF16(X) ((buf[(X-5)] << 8) | (buf[(X-4)]))
#define BUF32(X) ((buf[(X-5)] << 24) | (buf[(X-4)] << 16) | (buf[(X-3)] << 8) | (buf[(X-2)]))

#define GETPX(TID,LON) ((double)((LON) - tile[TID].lt_lon) / (tile[TID].rb_lon - tile[TID].lt_lon) * 256.0)
#define GETPY(TID,LAT) ((double)(tile[TID].lt_lat - (LAT)) / (tile[TID].lt_lat - tile[TID].rb_lat) * 256.0)
#define SETPV(TID,X,Y,V) (*(tile[TID].p + (Y) * 256 + (X)) = (V))

int grib_read(char *filename);
int tile_init();
int tile_clear();
int tile_search(int lat, int lon);
int tile_fill(unsigned char *rbuf, int lat0, int lon0, int dlat, int dlon, int grdw, int grdh);
int tile_save_pgm(time_t ft);
int tile_save(time_t ft);
int extract(unsigned char *src_buf, int src_len, int maxv, int nbits, unsigned char *dst_buf);
int fill_clut(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a);
int hexdump(unsigned char *buf, int len);

char mypath[1024];
char basedir[1024];
unsigned char clut_r[256],clut_g[256],clut_b[256],clut_a[256];
int verbose = 0;

int main(int argc, char **argv)
{
  char *p;
  strcpy(mypath,*argv);
  if ((p = rindex(mypath,'/')) != NULL) {
    *p = '\0';
  }

  /* initialize raw image buffer */
  tile_init();
  strcpy(basedir,"spl");

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

int grib_read(char *filename)
{
  FILE *fi, *fm;
  unsigned char buf[256000];
  unsigned char rawdata[256000];
  char msg[128];
  char fn[1024];
  int stat;
  struct tm reftm;
  time_t reftime;
  int slen,dlen,grdw,grdh,lat0,lon0,lat1,lon1,dlat,dlon;
  int pnum,type,ft,current_ft;
  int nbits,V,M;
  int ibuf,i;
  
  if ((fi = fopen(filename,"r")) == NULL) {
    perror(filename);
    exit(1);
  }

  fread(buf,1,16,fi);		/* section 0 */
  if (strncmp((const char *)buf,"GRIB",4) != 0) {
    fprintf(stderr,"not grib format\n");
    return -1;
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
  strftime(msg,sizeof(msg),"%Y-%m-%d.%T",&reftm);
  stat = BUF8(20);
  if (stat == 1){
    fprintf(stderr, "%s is test product. process stopped.\n", filename);
    exit(1);
  }

  current_ft = -999;
  while (fread(&slen,4,1,fi) > 0) {
    if (slen == 0x37373737) {	// section 8 = termination
      tile_save(reftime + current_ft * 60);	// flush final ft
      break;
    }

    slen = ntohl(slen);	/* section 3 */
    fread(buf,1,slen-4,fi);
    dlen = BUF32(7);
    grdw = BUF32(31);
    grdh = BUF32(35);
    if (verbose){printf("grid number %d x %d = %d\n",grdw,grdh,dlen);}
    lat0 = BUF32(47);
    lon0 = BUF32(51);
    lat1 = BUF32(56);
    lon1 = BUF32(60);
    dlon = BUF32(64);
    dlat = BUF32(68);
    if (verbose){printf("%d,%d - %d,%d / %d,%d\n",lat0,lon0,lat1,lon1,dlat,dlon);}

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 4 */
    fread(buf,1,slen-4,fi);
    pnum = BUF8(11);
    // parameter number : 203 = prec level
    type = BUF8(12);
    ft = BUF32(19);
    if (ft < 0) {
      ft = 0x80000000 - ft;
    }
    ft += 5;	// adjust for observation end time.
    if (ft != current_ft) {
      if (current_ft > -10) {
	tile_save(reftime + current_ft * 60);
      }
      if (verbose){printf("type = %d,%d, %d(%d)\n",pnum,type,ft,current_ft);}
      current_ft = ft;
      tile_clear();
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
    fill_clut(buf,clut_r,clut_g,clut_b,clut_a);
    // printf("V = %x (%d) ; M = %x (%d)\n", V,V,M,M);

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 6 */
    fread(buf,1,slen-4,fi);

    fread(&slen,4,1,fi);	slen = ntohl(slen);	/* section 7 */
    fread(buf,1,slen-4,fi);
    // if (verbose){hexdump(buf,slen-4);}

    if (pnum == 214) {      // parameter number : 214 = error info
      continue;
    }

    if (extract(&buf[1], slen-5, V, nbits, rawdata) != dlen) {
      fprintf(stderr, "Error: unmatch data length %d\n", dlen);
      exit(1);
    }
    tile_fill(rawdata, lat0 + dlat/2, lon0 - dlon/2, dlat , dlon, grdw, grdh);
  }
  fclose(fi);

  sprintf(fn,"%s/info.txt",basedir);
  if ((fm = fopen(fn,"w")) == NULL) {
    perror(fn);
    exit(1);
  }
  fprintf(fm,"reftm\t%ld\n",reftime);
  fprintf(fm,"reftstr\t%s\n",msg);
  fclose(fm);
  return 0;
}



/*
 * tile
 */

int tile_init()
{
  int i;
  for (i = 0; i < tile_count; i++) {
    tile[i].p = (unsigned char *)malloc(256 * 256);
    memset(tile[i].p, 0x00, 256 * 256);
  }
  return 0;
}

int tile_clear()
{
  int i;
  for (i = 0; i < tile_count; i++) {
    memset(tile[i].p, 0x00, 256 * 256);
  }
  return 0;
}

int tile_search(int lat, int lon)
{
  int i;
  for (i = 0; i < tile_count; i++) {
    if ((tile[i].lt_lat >= lat) && (lat > tile[i].rb_lat) &&
        (tile[i].lt_lon <= lon) && (lon < tile[i].rb_lon)) {
      return i;
    }
  }
  return -1;
}

int tile_pset(int tid, int px, int py, unsigned char val)
{
  *(tile[tid].p + py * 256 + px) = val;
  return 0;
}

int tile_fill(unsigned char *rbuf, int lat0, int lon0, int dlat, int dlon, int grdw, int grdh)
{
  int i,j,x,y,p = 0;
  int lat,lon;
  int tid,px,py,qx,qy,tid2,tid3,tid4;
  unsigned char c;
  int norishiro = 88;

  lat = lat0;
  lon = lon0;
  if ((tid = tile_search(lat, lon)) < 0){return -1;}
  for (y = 0; y < grdh; y++) {
    py = GETPY(tid, lat);
    if ((py < 0) || (py >= 256)) {
      if ((tid = tile_search(lat, lon)) < 0){return -1;}
      px = GETPX(tid, lon);
      py = GETPY(tid, lat);
    }

    lon = lon0;
    for (x = 0; x < grdw; x++) {
      px = GETPX(tid, lon);
      if ((px < 0) || (px >= 256)) {
        if ((tid = tile_search(lat, lon)) < 0){return -1;}
        px = GETPX(tid, lon);
	py = GETPY(tid, lat);
      }
      qx = GETPX(tid, lon + dlon + norishiro);
      qy = GETPY(tid, lat - dlat - norishiro);
      if (px == qx) {
	qx++;
      }
      if (py == qy) {
	qy++;
      }

      c = rbuf[p];

      if ((qx <= 256) && (qy <= 256)) {
        for (j = py; j < qy; j++) {
	  for (i = px; i < qx; i++) {
	    SETPV(tid,i,j,c);
	  }
	}
      }
      else {
	px = GETPX(tid, lon);
	py = GETPY(tid, lat);
	qx = GETPX(tid, lon + dlon + norishiro);
	qy = GETPY(tid, lat - dlat - norishiro);
	px = (px < 0) ? 0 : px;
	py = (py < 0) ? 0 : py;
	qx = (qx > 256) ? 256 : qx;
	qy = (qy > 256) ? 256 : qy;
        for (j = py; j < qy; j++) {
	  for (i = px; i < qx; i++) {
	    SETPV(tid,i,j,c);
	  }
	}

	if ((tid2 = tile_search(lat,        lon + dlon)) < 0){return -1;}
	if ((tid3 = tile_search(lat - dlat, lon       )) < 0){return -1;}
	if ((tid4 = tile_search(lat - dlat, lon + dlon)) < 0){return -1;}

	if (tid != tid2) {
	  px = GETPX(tid2, lon);
	  py = GETPY(tid2, lat);
	  qx = GETPX(tid2, lon + dlon + norishiro);
	  qy = GETPY(tid2, lat - dlat - norishiro);
	  px = (px < 0) ? 0 : px;
	  py = (py < 0) ? 0 : py;
	  qx = (qx > 256) ? 256 : qx;
	  qy = (qy > 256) ? 256 : qy;
	  for (j = py; j < qy; j++) {
	    for (i = px; i < qx; i++) {
	      SETPV(tid2,i,j,c);
	    }
	  }
	}
	if (tid2 != tid3) {
	  px = GETPX(tid3, lon);
	  py = GETPY(tid3, lat);
	  qx = GETPX(tid3, lon + dlon + norishiro);
	  qy = GETPY(tid3, lat - dlat - norishiro);
	  px = (px < 0) ? 0 : px;
	  py = (py < 0) ? 0 : py;
	  qx = (qx > 256) ? 256 : qx;
	  qy = (qy > 256) ? 256 : qy;
	  for (j = py; j < qy; j++) {
	    for (i = px; i < qx; i++) {
	      SETPV(tid3,i,j,c);
	    }
	  }
	}
	if ((tid2 != tid4) || (tid3 != tid4)) {
	  px = GETPX(tid4, lon);
	  py = GETPY(tid4, lat);
	  qx = GETPX(tid4, lon + dlon + norishiro);
	  qy = GETPY(tid4, lat - dlat - norishiro);
	  px = (px < 0) ? 0 : px;
	  py = (py < 0) ? 0 : py;
	  qx = (qx > 256) ? 256 : qx;
	  qy = (qy > 256) ? 256 : qy;
	  for (j = py; j < qy; j++) {
	    for (i = px; i < qx; i++) {
	      SETPV(tid4,i,j,c);
	    }
	  }
	}
      }

      lon += dlon;
      p++;
    }
    lat -= dlat;
  }
  return 0;
}

int tile_save_pgm(time_t ft)
{
  char outf[256];
  int i;
  FILE *fp;
  for (i = 0; i < tile_count; i++) {
    sprintf(outf,"spl/%s.pgm",tile[i].id);
    printf("%ld, %s\n",ft, outf);
    fp = fopen(outf,"w");
    fprintf(fp,"P5\n256 256\n255\n");
    fwrite(tile[i].p,1,256*256,fp);
    fclose(fp);
  }
  return 0;
}

int tile_save(time_t ft)
{
  struct tm *t;
  char outd[1024];
  char outf[1024];
  FILE *fp;
  png_structp     png_ptr;
  png_infop       info_ptr;
  png_bytep *rowp;
  int i,x,y,r;
  struct stat dstat;
  unsigned char *p, *q;
  png_color palette[256];
  png_byte trans[256];

  t = gmtime(&ft);
  strftime(outd,sizeof(outd),"%Y-%m-%d-%H-%M_Tiles",t);
  sprintf(outf,"%s/%s", basedir, outd);
  printf("output %s\n",outd);
  if (stat(outf,&dstat) < 0) {
    sprintf(outf,"mkdir -p %s/%s", basedir, outd);
    system(outf);
  }

  rowp = (png_bytepp)malloc(256 * sizeof(png_bytep));
  for (y = 0; y < 256; y++){
    rowp[y] = (unsigned char *)malloc(256 * 4);
  }

  for (i = 0; i < 256; i++) {
    palette[i].red   = clut_r[i];
    palette[i].green = clut_g[i];
    palette[i].blue  = clut_b[i];
    trans[i] = clut_a[i];
  }

  for (i = 0; i < tile_count; i++) {
    sprintf(outf,"%s/%s/%s.png", basedir, outd, tile[i].id);
    if (verbose){printf("saving tile %s\n",outf);}
    if ((fp = fopen(outf,"w")) == NULL) {
      perror(outf);
      return -1;
    }

    p = tile[i].p;
    for (y = 0; y < 256; y++){
      q = rowp[y];
      for (x = 0; x < 256; x++){
        *q++ = *p++;
      }
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
      fclose(fp);
      return -1;
    }
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, 256, 256, 8,
                 PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_PLTE(png_ptr, info_ptr, palette, 256);
    png_set_tRNS(png_ptr, info_ptr, trans, 256, NULL);
    png_set_rows(png_ptr, info_ptr, rowp);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
  }

  for (y = 0; y < 256; y++){
    free(rowp[y]);
  }
  free(rowp);

  sprintf(outf,"%s/grib-tile-post.sh %s/%s &", mypath, basedir, outd);
  system(outf);

  return 0;
}

/*
 * extract run-length
hmm...
 */
int extract(unsigned char *src_buf, int src_len, int maxv, int nbits, unsigned char *dst_buf)
{
  int lngu;
  int lastv;
  int pos = 0;
  int i;
  int lngv,lngn;

  lngu = (1 << nbits) - 1 - maxv;
  if (verbose) {
    printf("src_len = %d, maxv = %d, nbits = %d, lngu = %d\n", src_len, maxv, nbits, lngu);
  }

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

/*
 * clut : color lookup table
 */
int fill_clut(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a)
{
  int i, m, M, ibuf;
  unsigned int a = 256 * 0.55;
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
      clut_r[m] = clut_g[m] = clut_b[m] = 0xe0;
      clut_a[m] = a;
    }
    else if (ibuf < 400) {	// 1 .. 4 mm
      clut_r[m] = 0; clut_g[m] = 0x99; clut_b[m] = 255; clut_a[m] = a;
    }
    else if (ibuf < 1600) {	// 4 .. 16 mm
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0xcc; clut_a[m] = a;
    }
    else if (ibuf < 3200) {	// 16 .. 32 mm
      clut_r[m] = 255; clut_g[m] = 255; clut_b[m] = 0; clut_a[m] = a;
    }
    else if (ibuf < 30000) {	// 32 .. mm
      clut_r[m] = 255; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = a;
    }
    else {
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
    }
  }
  //m = 1;
  //clut_a[m] = 0;
  return 0;
}

int fill_clut_jma(unsigned char *buf, unsigned char *clut_r, unsigned char *clut_g, unsigned char *clut_b, unsigned char *clut_a)
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
      clut_r[m] = 180; clut_g[m] = 180; clut_b[m] = 180; clut_a[m] = 192;
    }
    else if (ibuf < 500) {	// 1 .. 5 mm
      clut_r[m] = 0; clut_g[m] = 216; clut_b[m] = 235; clut_a[m] = 192;
    }
    else if (ibuf < 1000) {	//  .. 10 mm
      clut_r[m] = 0; clut_g[m] = 160; clut_b[m] = 235; clut_a[m] = 192;
    }
    else if (ibuf < 2000) {	//  .. 20 mm
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 235; clut_a[m] = 192;
    }
    else if (ibuf < 3000) {	//  .. 30 mm
      clut_r[m] = 235; clut_g[m] = 235; clut_b[m] = 0; clut_a[m] = 192;
    }
    else if (ibuf < 5000) {	//  .. 50 mm
      clut_r[m] = 235; clut_g[m] = 192; clut_b[m] = 0; clut_a[m] = 192;
    }
    else if (ibuf < 8000) {	//  .. 80 mm
      clut_r[m] = 235; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 192;
    }
    else if (ibuf < 20000) {	// 32 .. mm
      clut_r[m] = 235; clut_g[m] = 0; clut_b[m] = 200; clut_a[m] = 192;
    }
    else {
      clut_r[m] = 0; clut_g[m] = 0; clut_b[m] = 0; clut_a[m] = 0;
    }
  }
  //m = 1;
  //clut_a[m] = 0;
  return 0;
}


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
