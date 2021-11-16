// TODO: for signed values, add 0x8000 first to make unsigned?
// Or shift sign bit (1s complement instead of 2s?)

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "rANS_static4x16.h"
#include "varint.h"

// 11.1s down to 5.9s and 253376580 up to 253583161.  Worth it.
//#define ED1 7
//#define ED2 2048
#define ED1 255
#define ED2 65536

#ifndef DECODE_ONLY
// Fast log approximation
double log_(double a) {
  union { double d; long long x; } u = { a };
  return (u.x - 4606921278410026770) * 1.539095918623324e-16;
}

// Byte-wise order-1 entropy calculation
double entropy8o1(unsigned char *data, int len) {
  int F[256][256] = {0}, T[256] = {0};
  double e = 0;
  int i, j;

  for (j = i = 0; i < len; i++) {
    F[j][data[i]]++;
    T[j]++;
    j=data[i];
    if ((i&ED1)==ED1) i+=ED2;
  }

  for (j = 0; j < 256; j++) {
    if (!T[j])
      continue;
    for (i = 0; i < 256; i++) {
      if (F[j][i]) {
	e += -log_((double)F[j][i]/T[j]) * F[j][i];
      }
    }
  }

  return e / log_(256);
}

// Striped order-1 entropy calculation
double entropy8o1b(unsigned char *data, int len, int stride) {
  int F[256][256] = {0}, T[256] = {0};
  double e = 0;
  int i, j;

  for (j = i = stride; i < len; i++) {
    F[data[i-stride]][data[i]]++;
    T[data[i-stride]]++;
    if ((i&ED1)==ED1) i+=ED2;
  }

  for (j = 0; j < 256; j++) {
    if (!T[j])
      continue;
    for (i = 0; i < 256; i++) {
      if (F[j][i]) {
	e += -log_((double)F[j][i]/T[j]) * F[j][i];
      }
    }
  }

  return e / log_(256);
}

// Compare each image to the previous image.
// Plus each row to previous row for 1st image
void delta2(char *buf, int len, int bits, int x, int y) {
  int i, j, xy = x*y;
  
  if (bits == 2) {
    uint16_t *buf_i = (uint16_t *)buf;
    len /= 2;

    // image by image from image 1 onwards
    for (i = len-xy; i > 0; i-=xy) {
      uint16_t *buf_i0 = &buf_i[i];
      uint16_t *buf_i1 = &buf_i[i-xy];
      for (j = 0; j < xy; j++)
	buf_i0[j] = buf_i0[j] - buf_i1[j];// + 0x8080;
    }

    // row by row for 1st image; edit: all images!
    //for (i = xy-x; i > 0; i-=x)
    for (i = len-x; i > 0; i-=x) {
      uint16_t *buf_i0 = &buf_i[i];
      uint16_t *buf_i1 = &buf_i[i-x];
      for (j = 0; j < x; j++)
	buf_i0[j] = buf_i0[j] - buf_i1[j] + 0x8080;
    }
  } else if (bits == 4) {
    uint32_t *buf_i = (uint32_t *)buf;

    len /= 4;
    for (i = len-xy; i > 0; i-=xy)
      for (j = 0; j < xy; j++)
	buf_i[i+j] = buf_i[i+j] - buf_i[i+j-xy];// + 0x80808080;

    // row by row
    //for (i = xy-x; i > 0; i-=x)
    for (i = len-x; i > 0; i-=x)
      for (j = 0; j < x; j++)
	buf_i[i+j] = buf_i[i+j] - buf_i[i+j-x] + 0x80808080;
  }
}
#endif // DECODE_ONLY

void undelta2(char *buf, int len, int bits, int x, int y) {
  int i, j, xy = x*y;

  fprintf(stderr, "len=%d b=%d x=%d y=%d\n", len, bits, x, y);
  
  if (bits == 2) {
    uint16_t *buf_i = (uint16_t *)buf;
    len /= 2;

    // row by row for 1st image; edit all images!
    //for (i = x; i < xy; i+=x)
    for (i = x; i < len; i+=x) {
      uint16_t *buf_i0 = &buf_i[i];
      uint16_t *buf_i1 = &buf_i[i-x];
      for (j = 0; j < x; j++)
	buf_i0[j] += buf_i1[j] - 0x8080;
    }

    // then image by image delta
    for (i = xy; i < len; i+=xy) {
      uint16_t *buf_i0 = &buf_i[i];
      uint16_t *buf_i1 = &buf_i[i-xy];
      for (j = 0; j < xy; j++)
	buf_i0[j] += buf_i1[j]; // -0x8080;
    }

  } else if (bits == 4) {
    uint32_t *buf_i = (uint32_t *)buf;
    len /= 4;

    //for (i = x; i < xy; i+=x)
    for (i = x; i < len; i+=x)
      for (j = 0; j < x; j++)
	buf_i[i+j] = buf_i[i+j] + buf_i[i+j-x] - 0x80808080;

    for (i = xy; i < len; i+=xy)
      for (j = 0; j < xy; j++)
	buf_i[i+j] = buf_i[i+j] + buf_i[i+j-xy];// - 0x80808080;
  }
}


#define CTX_SIZE 12
#define ZIG_ZAG

#ifndef DECODE_ONLY
#ifndef ZIG_ZAG
// TODO: get zig-zag decode done so we can use that instead
void delta_float(char *buf, int len, int bits, int x, int y) {
  unsigned int last[1<<CTX_SIZE] = {0};
  int i;

  // Delta vs previous value (in X)
  uint32_t *dat = (uint32_t *)buf;
  len /= 4;
  for (i = 0; i < len; i++) {
    uint32_t e = (dat[i]>>(32-CTX_SIZE));
    uint32_t m = dat[i];// & mask;
    uint32_t z = (m - last[e]); // & mask
    dat[i] = e | (z<<CTX_SIZE);
    last[e] = m;
  }
}
#else
// zig zag, but only 0.05% smaller on average?
void delta_float(char *buf, int len, int bits, int x, int y) {
  unsigned int last[1<<CTX_SIZE] = {0};
  int i;
  uint32_t *dat  = (uint32_t *)buf;

  // Delta vs previous value (in X)
  len /= 4;
  y = len/x;
  int iy, ix;
  for (iy = 0; iy < y; iy++) {
    for (ix = 0; ix < x; ix++) {
      i = iy*x+ix;
      uint32_t e = (dat[i]>>(32-CTX_SIZE));
      uint32_t m = dat[i];// & mask;
      uint32_t z = (m - last[e]); // & mask
      dat[i] = e | (z<<CTX_SIZE);
      last[e] = m;
    }
    if (++iy >= y)
      break;
    //for (ix = 0; ix < x; ix++) {
    for (ix = x-1; ix >= 0; ix--) {
      i = iy*x+ix;
      uint32_t e = (dat[i]>>(32-CTX_SIZE));
      uint32_t m = dat[i];// & mask;
      uint32_t z = (m - last[e]); // & mask
      dat[i] = e | (z<<CTX_SIZE);
      last[e] = m;
    }
  }
}
#endif // ZIG_ZAG

#endif // DECODE_ONLY

#ifndef ZIG_ZAG
// normal
void undelta_float(char *buf, int len, int bits, int x, int y) {
  unsigned int last[1<<CTX_SIZE] = {0};
  uint32_t *dat = (uint32_t *)buf;
  const uint32_t mask = (1<<(32-CTX_SIZE))-1;
  int i;

  len /= 4;
  for (i = 0; i < len; i++) {
    uint32_t e = dat[i] & ((1<<CTX_SIZE)-1);
    uint32_t m = dat[i]>>CTX_SIZE;
    m = (m + last[e]) & mask;
    last[e] = m;
    dat[i] = (e << (32-CTX_SIZE)) | m;
  }
}
#else
// zig zag
void undelta_float(char *buf, int len, int bits, int x, int y) {
  unsigned int last[1<<CTX_SIZE] = {0};
  uint32_t *dat = (uint32_t *)buf;
  const uint32_t mask = (1<<(32-CTX_SIZE))-1;
  int i;

  len /= 4;
  y = len/x;
  int iy, ix;

  for (iy = 0; iy < y; iy++) {
    for (ix = 0; ix < x; ix++) {
      i = iy*x+ix;
      uint32_t e = dat[i] & ((1<<CTX_SIZE)-1);
      uint32_t m = dat[i]>>CTX_SIZE;
      m = (m + last[e]) & mask;
      last[e] = m;
      dat[i] = (e << (32-CTX_SIZE)) | m;
    }
    if (++iy >= y)
      break;
    for (ix = x-1; ix >= 0; ix--) {
      i = iy*x+ix;
      uint32_t e = dat[i] & ((1<<CTX_SIZE)-1);
      uint32_t m = dat[i]>>CTX_SIZE;
      m = (m + last[e]) & mask;
      last[e] = m;
      dat[i] = (e << (32-CTX_SIZE)) | m;
    }
  }
}
#endif

#ifndef DECODE_ONLY
void encode(FILE *outfp, char *fn, unsigned char *buf, int len) {
  int known = 0;

  // Parse name: <NAME>_<X>x<Y>x<COUNT>x<BITS><FMT>
  char *cp = strrchr(fn, '_');
  if (!cp) goto unknown;

  int x = strtol(cp+1, &cp, 10);
  if (!cp) goto unknown;

  int y = strtol(cp+1, &cp, 10);
  if (!cp) goto unknown;

  int c = strtol(cp+1, &cp, 10);
  if (!cp) goto unknown;

  int b = strtol(cp+1, &cp, 10) / 8;
  if (!cp) goto unknown;

  char fmt = *cp;
  known = 1;

  if (fmt == 's') {
    int i;
    for (i = 0; i < len; i+=2)
      *(uint16_t *)&buf[i] = *(uint16_t *)&buf[i] + 0x8080;
  }

  unsigned char *buf2 = malloc(len);
  memcpy(buf2, buf, len);
  //    unsigned char *buf3 = malloc(len);
  //    memcpy(buf3, buf, len);
  unsigned char *cbuf = malloc(len*1.1);

  int i, j;
  unsigned char *out = buf;

  // buf  = orig.
  // bug2 = deltad version of buf
  if (fmt!='f') {
    delta2((char *)buf2, len, b, x, y);
    //delta3((char *)buf3, len, b, x, y);
  }
  else if (fmt == 'f')
    delta_float((char *)buf2, len, b, x, y);

  // Compute entropy of buf / buf2, with/without stripe (every 2/4 bytes).
  double e0 = entropy8o1b(buf,  len, b);
  double e2 = entropy8o1b(buf2, len, b);
  //    double e3 = entropy8o1b(buf3, len, b);
  //double E0 = e0 < e2 ? entropy8o1 (buf,  len) : e0+1;
  double E0 = e0 < e2 || fmt == 'f' ? entropy8o1 (buf,  len) : e0+1;
  fprintf(stderr, "Entropy before: %f %f %f\n", entropy8o1(buf,  len), e0, E0);
  fprintf(stderr, "Entropy delta2: %f %f b=%d\n", entropy8o1(buf2, len), e2, b);
  //    fprintf(stderr, "Entropy delta3: %f %f b=%d\n", entropy8o1(buf3, len), e3, b);

  // Pick encoding mode that works best
  int m;
  if (e0 < e2 || E0 < e2)
    m = 0, out = buf;
  else
    m = 2+(fmt=='f'), out = buf2;
  fprintf(stderr, "Choosing method m %d\n", m);

  // Header
  fputc(m + (fmt=='s')*4, outfp);
  fputc(b, outfp);
  fputc(x & 0xff, outfp); fputc(x >> 8, outfp);
  fputc(y & 0xff, outfp); fputc(y >> 8, outfp);
  fputc(len, outfp);
  fputc(len>>8, outfp);
  fputc(len>>16, outfp);
  fputc(len>>24, outfp);
  fflush(outfp);

 unknown:
  if (known == 0)
    // FIXME: write "unknown" header
    fprintf(stderr, "unknown\n");

  // Choose rANS "order" field, based on stats above
  int order = 9;
  switch(m) {
  case 0: order = E0 < e0 ? 1 : (b<<8) + 9; break;
  case 2: order = (b<<8) + 9; break;
  }
  order |= 4 + RANS_ORDER_STRIPE_NO0;

  // higher level compression:
  //order |= 128;// + 64;
  //order |= 128+64;
  //order |= 64;

  // Choose a sensible block size for optimal compression
  int target_size = 3000000;
  if (!(order & 8))
    target_size /= b;
  if (fmt == 'f')
    target_size *= (order & 8) ? 2 : 8; //243130230; 9.78, 9.84, 9.81
  while (x*y*b > target_size) {
    if (b%2 == 0 && b > 2)
      b /= 2;
    else if (x%2 == 0 && x > 4)
      x /= 2;
    else if (y%2 == 0 && y > 4)
      y /= 2;
    else
      break;
  }

  while (x*y*b*2 < target_size && x*y*b*2 < len)
    b *= 2;

  // Finally entropy encode
  while (len) {
    int len2;
    int ll = x*y*b > len ? len : x*y*b;
    len2 = rans_compress_bound_4x16(ll, order);
    char *comp = (char *)rans_compress_to_4x16(out, ll/*x*y*b*/, cbuf,
					       (unsigned int *)&len2, order);
    fprintf(stderr, "Compress %d to %d with order %x, remaining %d\n", x*y*b, len2, order, len-ll);

    unsigned char sz[4];
    sz[0] = len2;
    sz[1] = len2>>8;
    sz[2] = len2>>16;
    sz[3] = len2>>24;
    if (fwrite(sz, 1, 4, outfp) != 4) abort();
    if (fwrite(comp, 1, len2, outfp) != len2) abort();

    out += ll;
    len -= ll;
  }

  free(buf2);
  free(cbuf);
}
#endif // DECODE_ONLY

void decode(FILE *outfp, unsigned char *buf,  int len) {
  unsigned char *in = buf; len -= 10;
  int m = *in++;
  int b = *in++;
  int x = in[0] + in[1]*256; in+=2;
  int y = in[0] + in[1]*256; in+=2;
  unsigned int ulen = in[0] + (in[1]<<8) + (in[2]<<16) + (in[3]<<24);
  in += 4;
  fprintf(stderr, "m=%d b=%d x=%d y=%d, ulen=%d\n", m, b, x, y, ulen);
  char *out = malloc(ulen), *outp = out;
  if (!out)
    exit(1);

  char *ubuf = malloc(ulen);
    
  int ulen_tmp = ulen;
  while (ulen) {
    uint32_t clen = in[0] + (in[1]<<8) + (in[2]<<16) + (in[3]<<24), len2;
    in += 4;
    fprintf(stderr, "input %d\n", clen);
    var_get_u32(in+1, in+clen-1, &len2);
    unsigned char *uncomp = rans_uncompress_to_4x16(in, clen,
						    (unsigned char *)ubuf,
						    &len2);
    fprintf(stderr, "len2 = %d vs %d\n", len2, x*y*b);
    memcpy(outp, uncomp, len2); outp += len2;
    ulen -= len2;
    in += clen;
  }
  len = ulen_tmp;

  switch (m & 3) {
  case 0:
    break; // do nothing.
  case 2:
    undelta2(out, len, b, x, y);
    break;
  case 3:
    undelta_float(out, len, b, x, y);
    break;
  default:
    abort();
  }

  if (m & 4) {
    int i;
    for (i = 0; i < len; i+=2)
      *(uint16_t *)&out[i] =  *(uint16_t *)&out[i] - 0x8080;
  }

  if (fwrite(out, 1, len, outfp) != len)
    exit(1);

  free(out);
  free(ubuf);
}

int main(int argc, char **argv) {
  if (argc < 4)
    exit(1);

  // Load input file
  int enc = *argv[1] == 'e' ? 1 : 0;
  FILE *outfp = fopen(argv[3], "wb");
  if (!outfp) {
    perror(argv[3]);
    exit(1);
  }

  struct stat sb;
  if (stat(argv[2], &sb) < 0) {
    perror(argv[2]);
    exit(1);
  }
  int len = sb.st_size;

  unsigned char *buf = malloc(len);
  int fd = open(argv[2], O_RDONLY);
  if (fd < 0) {
    perror(argv[2]);
    exit(1);
  }
  if (read(fd, buf, len) != len)
    exit(1);
  close(fd);

  // do it
  if (enc)
#ifndef DECODE_ONLY
    encode(outfp, argv[2], buf, len);
#else
    exit(1);
#endif // DECODE_ONLY
  else
    decode(outfp, buf, len);

  return 0;
}
