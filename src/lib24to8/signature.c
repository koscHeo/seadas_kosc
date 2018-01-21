/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%        SSSSS  IIIII   GGGG  N   N   AAA   TTTTT  U   U  RRRR   EEEEE        %
%        SS       I    G      NN  N  A   A    T    U   U  R   R  E            %
%         SSS     I    G  GG  N N N  AAAAA    T    U   U  RRRR   EEE          %
%           SS    I    G   G  N  NN  A   A    T    U   U  R R    E            %
%        SSSSS  IIIII   GGG   N   N  A   A    T     UUU   R  R   EEEEE        %
%                                                                             %
%                                                                             %
%               Compute a Digital Signature for a Image Colormap              %
%                                                                             %
%                                                                             %
%                                                                             %
%                           Software Design                                   %
%                             John Cristy                                     %
%                            December 1992                                    %
%                                                                             %
%  Copyright 1994 E. I. du Pont de Nemours and Company                        %
%                                                                             %
%  Permission to use, copy, modify, distribute, and sell this software and    %
%  its documentation for any purpose is hereby granted without fee,           %
%  provided that the above Copyright notice appear in all copies and that     %
%  both that Copyright notice and this permission notice appear in            %
%  supporting documentation, and that the name of E. I. du Pont de Nemours    %
%  and Company not be used in advertising or publicity pertaining to          %
%  distribution of the software without specific, written prior               %
%  permission.  E. I. du Pont de Nemours and Company makes no representations %
%  about the suitability of this software for any purpose.  It is provided    %
%  "as is" without express or implied warranty.                               %
%                                                                             %
%  E. I. du Pont de Nemours and Company disclaims all warranties with regard  %
%  to this software, including all implied warranties of merchantability      %
%  and fitness, in no event shall E. I. du Pont de Nemours and Company be     %
%  liable for any special, indirect or consequential damages or any           %
%  damages whatsoever resulting from loss of use, data or profits, whether    %
%  in an action of contract, negligence or other tortuous action, arising     %
%  out of or in connection with the use or performance of this software.      %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Routine ColormapSignature computes a digital signature from the image
%  colormap.  This signature uniquely identifies the colormap and is convenient
%  for determining if the colormap of a sequence of images is identical when
%  animating.  The digital signature is from RSA Data Security MD5 Digest
%  Algorithm described in Internet draft [MD5], July 1992.
%
%
*/

/*
  Include declarations.
*/
#include "magick.h"
#include "image.h"

/*
  Typedef declarations.
*/
typedef struct _MessageDigest
{
  uint32_t 
    number_bits[2],
    accumulator[4];

  unsigned char
    message[64],
    digest[16];
} MessageDigest;

/*
  Forward declarations.
*/
static void
  TransformMessageDigest _Declare((MessageDigest *,uint32_t *)),
  UpdateMessageDigest _Declare((MessageDigest *,unsigned char *,uint32_t));

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C o m p u t e M e s s a g e D i g e s t                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function ComputeMessageDigest computes the message digest.
%
%  The format of the ComputeMessageDigest routine is:
%
%      ComputeMessageDigest(message_digest,input_message,message_length)
%
%  A description of each parameter follows:
%
%    o message_digest: The address of a structure of type MessageDigest.
%
%
*/
static void ComputeMessageDigest(message_digest)
MessageDigest
  *message_digest;
{
  int
    number_bytes;

  register unsigned char
    *p;

  register unsigned int
    i;

  unsigned char
    padding[64];

  uint32_t 
    message[16],
    padding_length;

  /*
    Save number of bits.
  */
  message[14]=message_digest->number_bits[0];
  message[15]=message_digest->number_bits[1];
  /*
    Compute number of bytes mod 64.
  */
  number_bytes=(int) ((message_digest->number_bits[0] >> 3) & 0x3F);
  /*
    Pad message to 56 mod 64.
  */
  padding_length=(number_bytes < 56) ? (56-number_bytes) : (120-number_bytes);
  padding[0]=0x80;
  for (i=1; i < padding_length; i++)
    padding[i]=(char) 0;
  UpdateMessageDigest(message_digest,padding,padding_length);
  /*
    Append length in bits and transform.
  */
  p=message_digest->message;
  for (i=0; i < 14; i++)
  {
    message[i]=(uint32_t) (*p++);
    message[i]|=((uint32_t) (*p++)) << 8;
    message[i]|=((uint32_t) (*p++)) << 16;
    message[i]|=((uint32_t) (*p++)) << 24;
  }
  TransformMessageDigest(message_digest,message);
  /*
    Store message in digest.
  */
  p=message_digest->digest;
  for (i=0; i < 4; i++)
  {
    *p++=(unsigned char) (message_digest->accumulator[i] & 0xff);
    *p++=(unsigned char) ((message_digest->accumulator[i] >> 8) & 0xff);
    *p++=(unsigned char) ((message_digest->accumulator[i] >> 16) & 0xff);
    *p++=(unsigned char) ((message_digest->accumulator[i] >> 24) & 0xff);
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   I n i t i a l i z e M e s s a g e D i g e s t                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function InitializeMessageDigest initializes the message digest structure.
%
%  The format of the InitializeMessageDigest routine is:
%
%      InitializeMessageDigest(message_digest)
%
%  A description of each parameter follows:
%
%    o message_digest: The address of a structure of type MessageDigest.
%
%
*/
static void InitializeMessageDigest(message_digest)
MessageDigest
  *message_digest;
{
  message_digest->number_bits[0]=(uint32_t) 0;
  message_digest->number_bits[1]=(uint32_t) 0;
  /*
    Load magic initialization constants.
  */
  message_digest->accumulator[0]=(uint32_t) 0x67452301;
  message_digest->accumulator[1]=(uint32_t) 0xefcdab89;
  message_digest->accumulator[2]=(uint32_t) 0x98badcfe;
  message_digest->accumulator[3]=(uint32_t) 0x10325476;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   T r a n s f o r m M e s s a g e D i g e s t                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function TransformMessageDigest updates the message digest.
%
%  The format of the TransformMessageDigest routine is:
%
%      TransformMessageDigest(message_digest,message)
%
%  A description of each parameter follows:
%
%    o message_digest: The address of a structure of type MessageDigest.
%
%
*/
static void TransformMessageDigest(message_digest,message)
MessageDigest
  *message_digest;

uint32_t 
  *message;
{
#define F(x,y,z)  (((x) & (y)) | ((~x) & (z)))
#define G(x,y,z)  (((x) & (z)) | ((y) & (~z)))
#define H(x,y,z)  ((x) ^ (y) ^ (z))
#define I(x,y,z)  ((y) ^ ((x) | (~z)))
#define RotateLeft(x,n)  (((x) << (n)) | (((x) & 0xffffffff) >> (32-(n))))

  static uint32_t 
    additive_constant[64]=  /* 4294967296*abs(sin(i)), i in radians */
    {
      0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf,
      0x4787c62a, 0xa8304613, 0xfd469501, 0x698098d8, 0x8b44f7af,
      0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e,
      0x49b40821, 0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
      0xd62f105d, 0x2441453,  0xd8a1e681, 0xe7d3fbc8, 0x21e1cde6,
      0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8,
      0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122,
      0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
      0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x4881d05,  0xd9d4d039,
      0xe6db99e5, 0x1fa27cf8, 0xc4ac5665, 0xf4292244, 0x432aff97,
      0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d,
      0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
      0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391,
    };

  register int
    i;

  register unsigned int
    j;

  register uint32_t 
    a,
    b,
    c,
    d,
    *p;

  /*
    Save accumulator.
  */
  a=message_digest->accumulator[0];
  b=message_digest->accumulator[1];
  c=message_digest->accumulator[2];
  d=message_digest->accumulator[3];
  /*
    a=b+((a+F(b,c,d)+X[k]+t) <<< s).
  */
  p=additive_constant;
  j=0;
  for (i=0; i < 4; i++)
  {
    a+=F(b,c,d)+message[j & 0x0f]+(*p++);
    a=RotateLeft(a,7)+b;
    j++;
    d+=F(a,b,c)+message[j & 0x0f]+(*p++);
    d=RotateLeft(d,12)+a;
    j++;
    c+=F(d,a,b)+message[j & 0x0f]+(*p++);
    c=RotateLeft(c,17)+d;
    j++;
    b+=F(c,d,a)+message[j & 0x0f]+(*p++);
    b=RotateLeft(b,22)+c;
    j++;
  }
  /*
    a=b+((a+G(b,c,d)+X[k]+t) <<< s).
  */
  j=1;
  for (i=0; i < 4; i++)
  {
    a+=G(b,c,d)+message[j & 0x0f]+(*p++);
    a=RotateLeft(a,5)+b;
    j+=5;
    d+=G(a,b,c)+message[j & 0x0f]+(*p++);
    d=RotateLeft(d,9)+a;
    j+=5;
    c+=G(d,a,b)+message[j & 0x0f]+(*p++);
    c=RotateLeft(c,14)+d;
    j+=5;
    b+=G(c,d,a)+message[j & 0x0f]+(*p++);
    b=RotateLeft(b,20)+c;
    j+=5;
  }
  /*
    a=b+((a+H(b,c,d)+X[k]+t) <<< s).
  */
  j=5;
  for (i=0; i < 4; i++)
  {
    a+=H(b,c,d)+message[j & 0x0f]+(*p++);
    a=RotateLeft(a,4)+b;
    j+=3;
    d+=H(a,b,c)+message[j & 0x0f]+(*p++);
    d=RotateLeft(d,11)+a;
    j+=3;
    c+=H(d,a,b)+message[j & 0x0f]+(*p++);
    c=RotateLeft(c,16)+d;
    j+=3;
    b+=H(c,d,a)+message[j & 0x0f]+(*p++);
    b=RotateLeft(b,23)+c;
    j+=3;
  }
  /*
    a=b+((a+I(b,c,d)+X[k]+t) <<< s).
  */
  j=0;
  for (i=0; i < 4; i++)
  {
    a+=I(b,c,d)+message[j & 0x0f]+(*p++);
    a=RotateLeft(a,6)+b;
    j+=7;
    d+=I(a,b,c)+message[j & 0x0f]+(*p++);
    d=RotateLeft(d,10)+a;
    j+=7;
    c+=I(d,a,b)+message[j & 0x0f]+(*p++);
    c=RotateLeft(c,15)+d;
    j+=7;
    b+=I(c,d,a)+message[j & 0x0f]+(*p++);
    b=RotateLeft(b,21)+c;
    j+=7;
  }
  /*
    Increment accumulator.
  */
  message_digest->accumulator[0]+=a;
  message_digest->accumulator[1]+=b;
  message_digest->accumulator[2]+=c;
  message_digest->accumulator[3]+=d;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   U p d a t e M e s s a g e D i g e s t                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function UpdateMessageDigest updates the message digest.
%
%  The format of the UpdateMessageDigest routine is:
%
%      UpdateMessageDigest(message_digest,input_message,message_length)
%
%  A description of each parameter follows:
%
%    o message_digest: The address of a structure of type MessageDigest.
%
%
*/
static void UpdateMessageDigest(message_digest,input_message,message_length)
MessageDigest
  *message_digest;

unsigned char
  *input_message;

uint32_t 
  message_length;
{
  int
    number_bytes;

  register unsigned char
    *p;

  register unsigned int
    i;

  uint32_t 
    message[16];

  /*
    Compute number of bytes mod 64.
  */
  number_bytes=(int) ((message_digest->number_bits[0] >> 3) & 0x3F);
  /*
    Update number of bits.
  */
  if (((message_digest->number_bits[0]+(message_length << 3)) & 0xffffffff) <
      message_digest->number_bits[0])
    message_digest->number_bits[1]++;
  message_digest->number_bits[0]+=message_length << 3;
  message_digest->number_bits[1]+=message_length >> 29;
  while (message_length--)
  {
    /*
      Add new character to message.
    */
    message_digest->message[number_bytes++]=(*input_message++);
    if (number_bytes == 0x40)
      {
        /*
          Transform message digest 64 bytes at a time.
        */
        p=message_digest->message;
        for (i=0; i < 16; i++)
        {
          message[i]=(uint32_t) (*p++);
          message[i]|=((uint32_t) (*p++)) << 8;
          message[i]|=((uint32_t) (*p++)) << 16;
          message[i]|=((uint32_t) (*p++)) << 24;
        }
        TransformMessageDigest(message_digest,message);
        number_bytes=0;
      }
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C o l o r m a p S i g n a t u r e                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Fucntion ColormapSignature computes a digital signature from the image
%  colormap.  This signature uniquely identifies the colormap and is convenient
%  for determining if the colormap of a sequence of images is identical when
%  animating.  The digital signature is from RSA Data Security MD5 Digest
%  Algorithm described in Internet draft [MD5], July 1992.
%
%  The format of the ColormapSignature routine is:
%
%      ColormapSignature(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type Image.
%
%
%
*/
void ColormapSignature(image)
Image
  *image;
{
  char
    *hex;

  MessageDigest
    message_digest;

  register int
    i;

  register unsigned char
    *p;

  unsigned char
    *colormap;

  if (image->class != PseudoClass)
    return;
  /*
    Allocate colormap.
  */
  colormap=(unsigned char *) malloc(3*image->colors*sizeof(unsigned char));
  if (colormap == (unsigned char *) NULL)
    {
      Warning("Unable to compute colormap signature",
        "Memory allocation failed");
      return;
    }
  p=colormap;
  for (i=0; i < image->colors; i++)
  {
    *p++=image->colormap[i].red;
    *p++=image->colormap[i].green;
    *p++=image->colormap[i].blue;
  }
  /*
    Compute program colormap signature.
  */
  InitializeMessageDigest(&message_digest);
  UpdateMessageDigest(&message_digest,colormap,
    (uint32_t) (image->colors*3));
  ComputeMessageDigest(&message_digest);
  (void) free((char *) colormap);
  /*
    Allocate memory for signature.
  */
  if (image->signature != (char *) NULL)
    (void) free((char *) image->signature);
  image->signature=(char *) malloc(33*sizeof(char));
  if (image->signature == (char *) NULL)
    {
      Warning("Unable to compute colormap signature",
        "Memory allocation failed");
      return;
    }
  /*
    Convert digital signature to a 32 character hex string.
  */
  p=(unsigned char *) image->signature;
  hex="0123456789abcdef";
  for (i=0; i < 16; i++)
  {
    *p++=hex[(message_digest.digest[i] >> 4) & 0xf];
    *p++=hex[message_digest.digest[i] & 0xf];
  }
  *p='\0';
}
