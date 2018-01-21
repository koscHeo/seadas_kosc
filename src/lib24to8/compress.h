/*
  Compress utility routines.
*/
extern unsigned int
  BMPDecodeImage _Declare((unsigned char *,unsigned char *,unsigned int,
    unsigned int,unsigned int)),
  BMPEncodeImage _Declare((unsigned char *,unsigned char *,unsigned int,
    unsigned int)),
  HuffmanDecodeImage _Declare((Image *)),
  HuffmanEncodeImage _Declare((Image *)),
  LZWDecodeImage _Declare((Image *)),
  LZWEncodeFilter _Declare((FILE *,unsigned char *,unsigned int)),
  LZWEncodeImage _Declare((Image *,unsigned int)),
  PackbitsEncodeImage _Declare((Image *,unsigned char *,unsigned char *)),
  QDecodeImage _Declare((unsigned char *,unsigned char *,unsigned int,
    unsigned int)),
  QEncodeImage _Declare((unsigned char *,unsigned char *,unsigned int,
    unsigned int)),
  RunlengthDecodeImage _Declare((Image *)),
  RunlengthEncodeImage _Declare((Image *)),
  SUNDecodeImage _Declare((unsigned char *,unsigned char *,unsigned int,
    unsigned int));
