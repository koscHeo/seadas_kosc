#include <xmlUtils.h>

#include <genutils.h>

#define XML_STRING_LENGTH 512

/* -------------------------------------------------------------------------- */
int xmlGetAttributeStrOpt(node_t* node, const char* name, char* buffer, int bufferSize) {
    node_t* child = roxml_get_attr(node, (char*)name, 0);
    if(child == NULL) {
        return 0;
    }
    int size;
    roxml_get_content(child, buffer, bufferSize, &size);
    if(size == bufferSize) {
        printf("%s Line %d: Buffer not big enough to hold contents of attribute \"%s\" of:\n",
                __FILE__, __LINE__, name);
        xmlPrintNodeTree(node);
        exit(1);
    }
    trimBlanks(buffer);
    return 1;
}

/* -------------------------------------------------------------------------- */
int xmlGetAttributeIntOpt(node_t* node, const char* name, int* val) {
    char str[XML_STRING_LENGTH];
    char* ptr;
    int result = xmlGetAttributeStrOpt(node, name, str, XML_STRING_LENGTH);
    if(result == 0)
        return 0;
    long tmp = strtol(str, &ptr, 0);
    if(str != ptr) {
        *val = tmp;
        return 1;
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
int xmlGetChildStrOpt(node_t* node, const char* name, char* buffer, int bufferSize) {
    node_t* child = roxml_get_chld(node, (char*)name, 0);
    if(child == NULL) {
        return 0;
    }
    int size;
    roxml_get_content(child, buffer, bufferSize, &size);
    if(size == bufferSize) {
        printf("%s Line %d: Buffer not big enough to hold contents of:\n",
                __FILE__, __LINE__);
        xmlPrintNodeTree(child);
        exit(1);
    }
    trimBlanks(buffer);
    return 1;
}

/* -------------------------------------------------------------------------- */
int xmlGetChildIntOpt(node_t* node, const char* name, int* val) {
    char str[XML_STRING_LENGTH];
    char* ptr;
    int result = xmlGetChildStrOpt(node, name, str, XML_STRING_LENGTH);
    if(result == 0)
        return 0;
    long tmp = strtol(str, &ptr, 0);
    if(str != ptr) {
        *val = tmp;
        return 1;
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
int xmlGetChildFloatOpt(node_t* node, const char* name, float* val) {
    char str[XML_STRING_LENGTH];
    char* ptr;
    int result = xmlGetChildStrOpt(node, name, str, XML_STRING_LENGTH);
    if(result == 0)
        return 0;
    float tmp = strtof(str, &ptr);
    if(str != ptr) {
        *val = tmp;
        return 1;
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
int xmlGetChildDoubleOpt(node_t* node, const char* name, double* val) {
    char str[XML_STRING_LENGTH];
    char* ptr;
    int result = xmlGetChildStrOpt(node, name, str, XML_STRING_LENGTH);
    if(result == 0)
        return 0;
    double tmp = strtod(str, &ptr);
    if(str != ptr) {
        *val = tmp;
        return 1;
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
void xmlGetAttributeStr(node_t* node, const char* name, char* buffer, int bufferSize) {
    node_t* child = roxml_get_attr(node, (char*)name, 0);
    if(child == NULL) {
        printf("%s Line %d: \"%s\" is not an attribute of:\n", __FILE__, __LINE__, name);
        xmlPrintNodeTree(node);
        exit(1);
    }
    int size;
    roxml_get_content(child, buffer, bufferSize, &size);
    if(size == bufferSize) {
        printf("%s Line %d: Buffer not big enough to hold contents of:\n",
                __FILE__, __LINE__);
        xmlPrintNodeTree(child);
        exit(1);
    }
    trimBlanks(buffer);
}

/* -------------------------------------------------------------------------- */
double xmlGetAttributeDouble(node_t* node, const char* name) {
    char str[XML_STRING_LENGTH];
    xmlGetAttributeStr(node, name, str, XML_STRING_LENGTH);
    char* ptr;
    double val;
    val = strtod(str, &ptr);
    if(ptr == str) {
        printf("%s Line %d: attribute \"%s\" is not a number:\n",
                __FILE__, __LINE__, name);
        xmlPrintNodeTree(node);
        exit(1);
    }
    return val;
}

/* -------------------------------------------------------------------------- */
void xmlGetChildStr(node_t* node, const char* name, char* buffer, int bufferSize) {
    node_t* child = roxml_get_chld(node, (char*)name, 0);
    if(child == NULL) {
        printf("%s Line %d: \"%s\" is not a child of:\n", __FILE__, __LINE__, name);
        xmlPrintNodeTree(node);
        exit(1);
    }
    int size;
    roxml_get_content(child, buffer, bufferSize, &size);
    if(size == bufferSize) {
        printf("%s Line %d: Buffer not big enough to hold contents of:\n",
                __FILE__, __LINE__);
        xmlPrintNodeTree(child);
        exit(1);
    }
    trimBlanks(buffer);
}

/* -------------------------------------------------------------------------- */
double xmlGetChildDouble(node_t* node, const char* name) {
    char str[XML_STRING_LENGTH];
    xmlGetChildStr(node, name, str, XML_STRING_LENGTH);
    char* ptr;
    double val;
    val = strtod(str, &ptr);
    if(ptr == str) {
        printf("%s Line %d: child \"%s\" is not a number:\n", __FILE__, __LINE__, name);
        xmlPrintNodeTree(node);
        exit(1);
    }
    return val;
}

/* -------------------------------------------------------------------------- */
void xmlPrintNodeTree(node_t* node) {
    node_t* parent = roxml_get_parent(node);
    if(parent != roxml_get_root(node))
        xmlPrintNodeTree(parent);

    char* str = roxml_get_name(node, NULL, 0);
    printf("<%s", str);
    roxml_release(str);

    int i;
    node_t* attr;
    int num = roxml_get_attr_nb(node);
    for(i=0; i<num; i++) {
        attr = roxml_get_attr(node, NULL, i);

        str = roxml_get_name(attr, NULL, 0);
        printf(" %s=", str);
        roxml_release(str);

        str = roxml_get_content(attr, NULL, 0, NULL);
        printf("\"%s\"", str);
        roxml_release(str);
    }
    printf(">\n");
}
