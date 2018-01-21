#include <stdint.h>

/*
  Image define declarations.
*/
#define AbsoluteValue(x)  ((x) < 0 ? -(x) : (x))
#define DegreesToRadians(x) ((x)*3.14159265358979323846/180.0)
#define Intensity(color)  (unsigned int)  \
  ((unsigned int) ((color).red*77+(color).green*150+(color).blue*29) >> 8)
#define MaxColormapSize  65535
#define MaxImageSize  (4096*4096)
#define MaxRGB  255
#define MaxRunlength  255
#define MaxTextLength  2048
#define RadiansToDegrees(x) ((x)*180/3.14159265358979323846)

/*
  Image Id's
*/
#define UndefinedId  0
#define ImageMagickId  1
/*
  Image classes:
*/
#define UndefinedClass  0
#define DirectClass  1
#define PseudoClass  2
/*
  Image colorspaces:
*/
#define UndefinedColorspace  0
#define RGBColorspace  1
#define GRAYColorspace 2
#define OHTAColorspace  3
#define XYZColorspace  4
#define YCbCrColorspace  5
#define YIQColorspace  6
#define YUVColorspace  7
/*
  Image compression algorithms:
*/
#define UndefinedCompression  0
#define NoCompression  1
#define RunlengthEncodedCompression  2
#define QEncodedCompression  3
/*
  Image interlace:
*/
#define UndefinedInterlace  0
#define NoneInterlace  1
#define LineInterlace  2
#define PlaneInterlace  3
/*
  Image compositing operations:
*/
#define UndefinedCompositeOp  0
#define OverCompositeOp  1
#define InCompositeOp  2
#define OutCompositeOp  3
#define AtopCompositeOp  4
#define XorCompositeOp  5
#define PlusCompositeOp  6
#define MinusCompositeOp  7
#define AddCompositeOp  8
#define SubtractCompositeOp  9
#define DifferenceCompositeOp  10
#define ReplaceCompositeOp  11
/*
  Page geometries:
*/
#define PSDensityGeometry  "72x72"
#define PSPageGeometry  "612x792+18+94"
#define TextPageGeometry  "612x792+36+36"

/*
  Typedef declarations for the Display program.
*/
typedef struct _ColorPacket
{
  unsigned char
    red,
    green,
    blue,
    flags;

  unsigned short
    index;
} ColorPacket;

typedef struct _ImageInfo
{
  char
    filename[MaxTextLength],
    magick[12];

  unsigned int
    assert;

  char
    *server_name,
    *font,
    *geometry,
    *density,
    *page;

  unsigned int
    interlace,
    monochrome,
    quality,
    verbose;

  char
    *undercolor;
} ImageInfo;

typedef struct _RectangleInfo
{
  unsigned int
    width,
    height;

  int
    x,
    y;
} RectangleInfo;

typedef struct _RunlengthPacket
{
  unsigned char
    red,
    green,
    blue,
    length;

  unsigned short
    index;
} RunlengthPacket;

typedef struct _Image
{
  FILE
    *file;

  int
    status;

  char
    filename[MaxTextLength];

  int
    pipe;

  char
    magick[12],
    *comments,
    *label;

  unsigned int
    id,
    class,
    colorspace,
    matte,
    compression,
    columns,
    rows,
    colors,
    scene;

  char
    *montage,
    *directory;

  ColorPacket
    *colormap;

  char
    *signature;

  RunlengthPacket
    *pixels,
    *packet;

  uint32_t
     packets;

  unsigned int
    runlength,
    packet_size;

  unsigned char
    *packed_pixels;

  unsigned int
    orphan;

  struct _Image
    *previous,
    *next;
} Image;

/*
  Image utilities routines.
*/
extern void
  CommentImage _Declare((Image *,char *)),
  Error _Declare((char *,char *)),
  LabelImage _Declare((Image *,char *)),
  Warning _Declare((char *,char *));

extern Image
  *AllocateImage _Declare((ImageInfo *)),
  *BorderImage _Declare((Image *,RectangleInfo *,ColorPacket *)),
  *BlurImage _Declare((Image *)),
  *ClipImage _Declare((Image *,RectangleInfo *)),
  *CopyImage _Declare((Image *,unsigned int,unsigned int,unsigned int)),
  *CutImage _Declare((Image *,RectangleInfo *)),
  *DespeckleImage _Declare((Image *)),
  *EdgeImage _Declare((Image *)),
  *EnhanceImage _Declare((Image *)),
  *FlipImage _Declare((Image *)),
  *FlopImage _Declare((Image *)),
  *FrameImage _Declare((Image *,RectangleInfo *,unsigned int,ColorPacket *,
    ColorPacket *,ColorPacket *)),
  *NoisyImage _Declare((Image *)),
  *ReadImage _Declare((ImageInfo *)),
  *RollImage _Declare((Image *,int,int)),
  *RotateImage _Declare((Image *,double,ColorPacket *,unsigned int)),
  *SampleImage _Declare((Image *,unsigned int,unsigned int)),
  *ScaleImage _Declare((Image *,unsigned int,unsigned int)),
  *SharpenImage _Declare((Image *)),
  *ShearImage _Declare((Image *,double,double,ColorPacket *,unsigned int)),
  *StereoImage _Declare((Image *,Image *));

extern int
  ReadDataBlock _Declare((char *,FILE *));

extern unsigned int
  IsGrayImage _Declare((Image *)),
  NumberColors _Declare((Image *,FILE *)),
  ReadData _Declare((char *,int,int,FILE *)),
  UncompressImage _Declare((Image *)),
  WriteImage _Declare((ImageInfo *,Image *));

extern void
  CloseImage _Declare((Image *)),
  ColormapSignature _Declare((Image *)),
  CompositeImage _Declare((Image *,unsigned int,Image *,int,int)),
  CompressColormap _Declare((Image *)),
  CompressImage _Declare((Image *)),
  DestroyImage _Declare((Image *)),
  DestroyImages _Declare((Image *)),
  EqualizeImage _Declare((Image *)),
  GammaImage _Declare((Image *,char *)),
  GetImageInfo _Declare((char *,ImageInfo *)),
  NegateImage _Declare((Image *)),
  NormalizeImage _Declare((Image *)),
  OpenImage _Declare((Image *,char *)),
  ParseImageGeometry _Declare((char *,unsigned int *,unsigned int *)),
  QuantizationError _Declare((Image *,unsigned int *,double *,double *)),
  QuantizeImage _Declare((Image *,unsigned int,unsigned int,unsigned int,
    unsigned int,unsigned int)),
  QuantizeImages _Declare((Image **,unsigned int,Image *,unsigned int,
    unsigned int,unsigned int,unsigned int,unsigned int)),
  RGBTransformImage _Declare((Image *,unsigned int)),
  SetErrorHandler _Declare((ErrorHandler)),
  SetImageMagick _Declare((ImageInfo *)),
  SetWarningHandler _Declare((ErrorHandler)),
  SortColormapByIntensity _Declare((Image *)),
  SyncImage _Declare((Image *)),
  TransformImage _Declare((Image **,char *,char *)),
  TransformRGBImage _Declare((Image *,unsigned int));
