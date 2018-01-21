#ifndef XML_UTILS_H
#define XML_UTILS_H

#include <roxml.h>

#ifdef __cplusplus
extern "C" {
#endif

  int xmlGetAttributeStrOpt(node_t* node, const char* name, char* buffer, int bufferSize);
  int xmlGetAttributeIntOpt(node_t* node, const char* name, int* val);
  int xmlGetChildStrOpt(node_t* node, const char* name, char* buffer, int bufferSize);
  int xmlGetChildIntOpt(node_t* node, const char* name, int* val);
  int xmlGetChildFloatOpt(node_t* node, const char* name, float* val);
  int xmlGetChildDoubleOpt(node_t* node, const char* name, double* val);
  void xmlGetAttributeStr(node_t* node, const char* name, char* buffer, int bufferSize);
  double xmlGetAttributeDouble(node_t* node, const char* name);
  void xmlGetChildStr(node_t* node, const char* name, char* buffer, int bufferSize);
  double xmlGetChildDouble(node_t* node, const char* name);
  void xmlPrintNodeTree(node_t* node);

#ifdef __cplusplus
}
#endif

#endif
