#include <productInfo.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <genutils.h>
#include <roxml.h>
#include <xmlUtils.h>
#include <sensorInfo.h>

#define XML_STRING_SIZE 512

/*
 * global variables
 */
static char productXMLFileName[FILENAME_MAX];
static node_t* rootNode = NULL;
static node_t* productsNode;
static node_t* productNode;
static node_t* algorithmNode;

static productInfo_t* productInfo = NULL;
static productInfo_t* algorithmInfo = NULL;

/* 
 * dup string if src is not NULL.
 */
static char* duplicateString(const char* src) {
    if (src == NULL)
        return NULL;
    else
        return strdup(src);
}

/* 
 * copy string to dest if src is not equal to the defaultStr
 *
 * @return 1 if string copied, 0 if src equal to default
 */
//static int copyStringDefault(char** dest, const char* src, const char* defaultStr) {
//    if (strcmp(src, defaultStr)) {
//        if (*dest != NULL)
//            free(*dest);
//        if (src == NULL)
//            *dest = NULL;
//        else
//            *dest = strdup(src);
//        return 1;
//    } else
//        return 0;
//}

/* 
 * copy int to dest if src is not equal to the defaultInt
 */
//static void copyIntDefault(int* dest, int src, int defaultInt) {
//    if (src != defaultInt)
//        *dest = src;
//}

/* 
 * copy double to dest if src is not equal to the defaultInt
 */
//static void copyDoubleDefault(double* dest, double src, double defaultDouble) {
//    if (src != defaultDouble)
//        *dest = src;
//}

/**
 * set the product structure to defaults freeing memory from the current values
 *
 * @param info destination product structure
 */
void clearProductInfo(productInfo_t* info) {
    if(info->description)
      free(info->description);
    info->description = duplicateString(PRODUCT_DEFAULT_description);
    if(info->units)
      free(info->units);
    info->units = duplicateString(PRODUCT_DEFAULT_units);
    if(info->palette)
      free(info->palette);
    info->palette = duplicateString(PRODUCT_DEFAULT_palette);
    if(info->paramDesignator)
      free(info->paramDesignator);
    info->paramDesignator = duplicateString(PRODUCT_DEFAULT_paramDesignator);
    info->paramWaveMin = PRODUCT_DEFAULT_paramWaveMin;
    info->paramWaveMax = PRODUCT_DEFAULT_paramWaveMax;
    if(info->standardName)
      free(info->standardName);
    info->standardName = duplicateString(PRODUCT_DEFAULT_standardName);
    if(info->category)
      free(info->category);
    info->category = duplicateString(PRODUCT_DEFAULT_category);
    if(info->dataType)
      free(info->dataType);
    info->dataType = duplicateString(PRODUCT_DEFAULT_dataType);
    if(info->prefix)
      free(info->prefix);
    info->prefix = duplicateString(PRODUCT_DEFAULT_prefix);
    if(info->suffix)
      free(info->suffix);
    info->suffix = duplicateString(PRODUCT_DEFAULT_suffix);
    if(info->algorithmName)
      free(info->algorithmName);
    info->algorithmName = duplicateString(PRODUCT_DEFAULT_algorithmName);
    if(info->productName)
      free(info->productName);
    info->productName = duplicateString(PRODUCT_DEFAULT_productName);
    info->cat_ix = PRODUCT_DEFAULT_cat_ix;
    info->prod_ix = PRODUCT_DEFAULT_prod_ix;
    info->rank = PRODUCT_DEFAULT_rank;
    info->fillValue = PRODUCT_DEFAULT_fillValue;
    info->validMin = PRODUCT_DEFAULT_validMin;
    info->validMax = PRODUCT_DEFAULT_validMax;
    if(info->displayScale)
      free(info->displayScale);
    info->displayScale = duplicateString(PRODUCT_DEFAULT_displayScale);
    info->displayMin = PRODUCT_DEFAULT_displayMin;
    info->displayMax = PRODUCT_DEFAULT_displayMax;
    info->scaleFactor = PRODUCT_DEFAULT_scaleFactor;
    info->addOffset = PRODUCT_DEFAULT_addOffset;
    if(info->reference)
        free(info->reference);
    info->reference = duplicateString(PRODUCT_DEFAULT_reference);
    if(info->comment)
            free(info->comment);
    info->comment = duplicateString(PRODUCT_DEFAULT_comment);
}

/**
 * set the product structure to defaults ignoring the current values
 *
 * @param info destination product structure
 */
void initProductInfo(productInfo_t* info) {
    bzero(info, sizeof(productInfo_t));
    clearProductInfo(info);
}

/**
 * allocate memory for the product into structure and init to defaults.
 * structure should be freed using freeProductInfo()
 *
 * @return pointer to newly allocated structure
 */
productInfo_t* allocateProductInfo() {
    productInfo_t* info = (productInfo_t*) allocateMemory(sizeof(productInfo_t),
            "productInfo");
    clearProductInfo(info);
    return info;
}

/**
 * free all the internal memory and the productInfo structure memory.
 *
 * @param info pointer to the product structure
 */
void freeProductInfo(productInfo_t* info) {
    if(info->description)
      free(info->description);
    if(info->units)
      free(info->units);
    if(info->palette)
      free(info->palette);
    if(info->paramDesignator)
      free(info->paramDesignator);
    if(info->standardName)
      free(info->standardName);
    if(info->category)
      free(info->category);
    if(info->dataType)
      free(info->dataType);
    if(info->prefix)
      free(info->prefix);
    if(info->suffix)
      free(info->suffix);
    if(info->algorithmName)
      free(info->algorithmName);
    if(info->productName)
      free(info->productName);
    if(info->displayScale)
      free(info->displayScale);
    if(info->reference)
        free(info->reference);
    if(info->comment)
        free(info->comment);
    free(info);
}

/**
 * copy product info structure, just header data
 *
 * @param dest destination product structure
 * @param src source product structure
 */
void copyProductInfoHeader(productInfo_t* dest, const productInfo_t* src) {
    if(dest->productName)
      free(dest->productName);
    dest->productName = duplicateString(src->productName);
    if(dest->paramDesignator)
      free(dest->paramDesignator);
    dest->paramDesignator = duplicateString(src->paramDesignator);
    dest->paramWaveMin = src->paramWaveMin;
    dest->paramWaveMax = src->paramWaveMax;
}

/**
 * copy product info structure
 *
 * @param dest destination product structure
 * @param src source product structure
 */
void copyProductInfo(productInfo_t* dest, const productInfo_t* src) {
    if(dest->description && strlen(dest->description) > 0 )
      free(dest->description);
    dest->description = duplicateString(src->description);
    if(dest->units)
      free(dest->units);
    dest->units = duplicateString(src->units);
    if(dest->palette)
      free(dest->palette);
    dest->palette = duplicateString(src->palette);
    if(dest->paramDesignator)
      free(dest->paramDesignator);
    dest->paramDesignator = duplicateString(src->paramDesignator);
    dest->paramWaveMin = src->paramWaveMin;
    dest->paramWaveMax = src->paramWaveMax;
    if(dest->standardName)
      free(dest->standardName);
    dest->standardName = duplicateString(src->standardName);
    if(dest->category)
      free(dest->category);
    dest->category = duplicateString(src->category);
    if(dest->dataType)
      free(dest->dataType);
    dest->dataType = duplicateString(src->dataType);
    if(dest->prefix)
      free(dest->prefix);
    dest->prefix = duplicateString(src->prefix);
    if(dest->suffix)
      free(dest->suffix);
    dest->suffix = duplicateString(src->suffix);
    if(dest->algorithmName)
      free(dest->algorithmName);
    dest->algorithmName = duplicateString(src->algorithmName);
    if(dest->productName)
      free(dest->productName);
    dest->productName = duplicateString(src->productName);
    dest->cat_ix = src->cat_ix;
    dest->prod_ix = src->prod_ix;
    dest->rank = src->rank;
    dest->fillValue = src->fillValue;
    dest->validMin = src->validMin;
    dest->validMax = src->validMax;
    if(dest->displayScale)
      free(dest->displayScale);
    dest->displayScale = duplicateString(src->displayScale);
    dest->displayMin = src->displayMin;
    dest->displayMax = src->displayMax;
    dest->scaleFactor = src->scaleFactor;
    dest->addOffset = src->addOffset;
    if(dest->reference)
      free(dest->reference);
    dest->reference = duplicateString(src->reference);
    if(dest->comment)
        free(dest->comment);
      dest->comment = duplicateString(src->comment);
}

/**
 * add non default fields from src into dest
 *
 * @param dest destination product structure
 * @param src source product structure
 */
//void addProductInfo(productInfo_t* dest, const productInfo_t* src) {
//    copyStringDefault(&dest->name, src->name, PRODUCT_DEFAULT_name);
//    copyStringDefault(&dest->description, src->description,
//            PRODUCT_DEFAULT_description);
//    copyStringDefault(&dest->units, src->units, PRODUCT_DEFAULT_units);
//    copyStringDefault(&dest->palette, src->palette, PRODUCT_DEFAULT_palette);
//    if(copyStringDefault(&dest->paramDesignator, src->paramDesignator,
//            PRODUCT_DEFAULT_paramDesignator)) {
//        dest->paramWaveMin = src->paramWaveMin;
//        dest->paramWaveMax = src->paramWaveMax;
//    }
//    copyStringDefault(&dest->standardName, src->standardName,
//            PRODUCT_DEFAULT_standardName);
//    copyStringDefault(&dest->category, src->category, PRODUCT_DEFAULT_category);
//    copyStringDefault(&dest->dataType, src->dataType, PRODUCT_DEFAULT_dataType);
//    copyStringDefault(&dest->prefix, src->prefix, PRODUCT_DEFAULT_prefix);
//    copyStringDefault(&dest->suffix, src->suffix, PRODUCT_DEFAULT_suffix);
//    copyStringDefault(&dest->algorithmName, src->algorithmName,
//            PRODUCT_DEFAULT_algorithmName);
//    copyStringDefault(&dest->productName, src->productName,
//            PRODUCT_DEFAULT_productName);
//    copyIntDefault(&dest->cat_ix, src->cat_ix, PRODUCT_DEFAULT_cat_ix);
//    copyIntDefault(&dest->prod_ix, src->prod_ix, PRODUCT_DEFAULT_prod_ix);
//    copyIntDefault(&dest->rank, src->rank, PRODUCT_DEFAULT_rank);
//    copyDoubleDefault(&dest->fillValue, src->fillValue,
//            PRODUCT_DEFAULT_fillValue);
//    copyDoubleDefault(&dest->validMin, src->validMin, PRODUCT_DEFAULT_validMin);
//    copyDoubleDefault(&dest->validMax, src->validMax, PRODUCT_DEFAULT_validMax);
//    copyStringDefault(&dest->displayScale, src->displayScale,
//            PRODUCT_DEFAULT_displayScale);
//    copyDoubleDefault(&dest->displayMin, src->displayMin,
//            PRODUCT_DEFAULT_displayMin);
//    copyDoubleDefault(&dest->displayMax, src->displayMax,
//            PRODUCT_DEFAULT_displayMax);
//    copyDoubleDefault(&dest->scaleFactor, src->scaleFactor,
//            PRODUCT_DEFAULT_scaleFactor);
//    copyDoubleDefault(&dest->addOffset, src->addOffset,
//            PRODUCT_DEFAULT_addOffset);
//}


/**
 * load XML file into local XML structure.
 */
void initXmlFile() {

    // bail if root node is already set
    if(rootNode)
        return;

    char *dataRoot;
    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }
    strcpy(productXMLFileName, dataRoot);
    strcat(productXMLFileName,"/common/product.xml");

    rootNode = roxml_load_doc(productXMLFileName);
    if(rootNode == NULL) {
        printf("%s Line %d: could not open product XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        exit(1);
    }
    productsNode = roxml_get_chld(rootNode, "products", 0);
    if(productsNode == NULL) {
        printf("%s Line %d: could not find products tag in XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        exit(1);
    }

    productInfo = allocateProductInfo();
    algorithmInfo = allocateProductInfo();

}

/**
 * check the product node to make sure the name == "product"
 */
void checkProductNode() {
    char tmpStr[XML_STRING_SIZE];

    roxml_get_name(productNode, tmpStr, XML_STRING_SIZE);
    if(strcmp(tmpStr, "product") != 0) {
        printf("%s Line %d: \"product\" node expected, found \"%s\" in file %s\n",
                __FILE__, __LINE__, tmpStr, productXMLFileName);
        exit(1);
    }
}

/**
 * check the algorithm node to make sure the name == "algorithm"
 */
void checkAlgorithmNode() {
    char tmpStr[XML_STRING_SIZE];

    roxml_get_name(algorithmNode, tmpStr, XML_STRING_SIZE);
    if(strcmp(tmpStr, "algorithm") != 0) {
        printf("%s Line %d: \"algorithm\" node expected, found \"%s\" in file %s\n",
                __FILE__, __LINE__, tmpStr, productXMLFileName);
        exit(1);
    }
}

/**
 * read the paramDesignator information out of the XML file
 */
void readParamDesignator(productInfo_t* info, node_t* node) {
    char tmpStr[XML_STRING_SIZE];
    node_t* tmpNode;

    tmpNode = roxml_get_chld(node, "paramDesignator", 0);
    while(tmpNode) {

        // grab the text
        int size;
        roxml_get_content(tmpNode, tmpStr, XML_STRING_SIZE, &size);
        trimBlanks(tmpStr);

        if(strcmp(tmpStr, "none") == 0 ||
                strcmp(tmpStr, "band") == 0  ||
                strcmp(tmpStr, "int") == 0) {
            free(info->paramDesignator);
            info->paramDesignator = strdup(tmpStr);
            info->paramWaveMin = PRODUCT_DEFAULT_paramWaveMin;
            info->paramWaveMax = PRODUCT_DEFAULT_paramWaveMax;
            return;
        } else if(strcmp(tmpStr, "uv") == 0) {
            if(strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if(info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 100)
                info->paramWaveMin = 100;
            if(info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 400)
                info->paramWaveMax = 400;
        } else if(strcmp(tmpStr, "visible") == 0) {
            if(strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if(info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 400)
                info->paramWaveMin = 400;
            if(info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 725)
                info->paramWaveMax = 725;
        } else if(strcmp(tmpStr, "nir") == 0) {
            if(strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if(info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 725)
                info->paramWaveMin = 725;
            if(info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 1400)
                info->paramWaveMax = 1400;
        } else if(strcmp(tmpStr, "swir") == 0) {
            if(strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if(info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 1400)
                info->paramWaveMin = 1400;
            if(info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 3000)
                info->paramWaveMax = 3000;
        } else if(strcmp(tmpStr, "emissive") == 0) {
            if(strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if(info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 3000)
                info->paramWaveMin = 3000;
            if(info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 15000)
                info->paramWaveMax = 15000;
        } else {
            printf("%s Line %d: \"paramDesignator\" node has illegal value \"%s\" in file %s\n",
                    __FILE__, __LINE__, tmpStr, productXMLFileName);
            exit(1);
        }

        // grab next node
        tmpNode = roxml_get_next_sibling(tmpNode);

        // make sure the next node is another paramDesignator
        if(tmpNode) {
            roxml_get_name(tmpNode, tmpStr, XML_STRING_SIZE);
            if(strcmp(tmpStr, "paramDesignator") != 0)
                tmpNode = NULL;
        }
    }
}

/**
 * read the range node and write the data into the productInfo structure
 *
 * @param productInfo product structure to write the range data
 * @param rangeNode XML node to read the info from
 */
void readSingleRange(productInfo_t* productInfo, node_t* rangeNode) {
    xmlGetChildDoubleOpt(rangeNode, "validMin", &productInfo->validMin);
    xmlGetChildDoubleOpt(rangeNode, "validMax", &productInfo->validMax);
    xmlGetChildDoubleOpt(rangeNode, "displayMin", &productInfo->displayMin);
    xmlGetChildDoubleOpt(rangeNode, "displayMax", &productInfo->displayMax);
    xmlGetChildDoubleOpt(rangeNode, "scaleFactor", &productInfo->scaleFactor);
    xmlGetChildDoubleOpt(rangeNode, "addOffset", &productInfo->addOffset);
}

/**
 * read the range structures given the paramVal
 *
 * @param productInfo product structure to write the range data
 * @param productNode XML node to read the info from
 * @param paramVal integer to match to the range's min and max
 */
void readRange(productInfo_t* productInfo, node_t* productNode, int paramVal) {
    char tmpStr[XML_STRING_SIZE];
    int min, max;
    node_t* rangeNode;

    // set the default
    rangeNode = roxml_get_chld(productNode, "range", 0);
    do {

        // make sure it is a range node
        roxml_get_name(rangeNode, tmpStr, XML_STRING_SIZE);
        if(strcmp(tmpStr, "range") != 0)
            break;

        // make sure min and max are not defined
        if(xmlGetAttributeIntOpt(rangeNode, "min", &min))
            continue;
        if(xmlGetAttributeIntOpt(rangeNode, "max", &max))
            continue;

        readSingleRange(productInfo, rangeNode);

    } while((rangeNode = roxml_get_next_sibling(rangeNode)) != NULL);

    // now search for the matching range
    rangeNode = roxml_get_chld(productNode, "range", 0);
    do {

        // make sure it is a range node
        roxml_get_name(rangeNode, tmpStr, XML_STRING_SIZE);
        if(strcmp(tmpStr, "range") != 0)
            break;

        // make sure min and max are defined
        if(xmlGetAttributeIntOpt(rangeNode, "min", &min))
            if(xmlGetAttributeIntOpt(rangeNode, "max", &max))
                if(min<=paramVal && paramVal<=max)
                    readSingleRange(productInfo, rangeNode);

    } while((rangeNode = roxml_get_next_sibling(rangeNode)) != NULL);

}



/**
 * read product header from the XML productNode into productInfo
 */
void readProductHeader() {
    char tmpStr[XML_STRING_SIZE];

    // clear and set productName
    xmlGetAttributeStr(productNode, "name", tmpStr, XML_STRING_SIZE);
    if(productInfo->productName)
        free(productInfo->productName);
    productInfo->productName = strdup(tmpStr);

    // clear and set paramDesignator
    if(productInfo->paramDesignator)
      free(productInfo->paramDesignator);
    productInfo->paramDesignator = duplicateString(PRODUCT_DEFAULT_paramDesignator);
    readParamDesignator(productInfo, productNode);
}

/**
 * read product data from the XML productNode into productInfo
 */
void readProduct(int paramVal) {
    char tmpStr[XML_STRING_SIZE];

    // clear product structure to defaults
    clearProductInfo(productInfo);

    xmlGetAttributeStr(productNode, "name", tmpStr, XML_STRING_SIZE);
    if(productInfo->productName)
        free(productInfo->productName);
    productInfo->productName = strdup(tmpStr);

    if(xmlGetChildStrOpt(productNode, "standardName", tmpStr, XML_STRING_SIZE)) {
        if(productInfo->standardName)
            free(productInfo->standardName);
        productInfo->standardName = strdup(tmpStr);
    }

    xmlGetChildStr(productNode, "units", tmpStr, XML_STRING_SIZE);
    if(productInfo->units)
        free(productInfo->units);
    productInfo->units = strdup(tmpStr);

    if(xmlGetChildStrOpt(productNode, "palette", tmpStr, XML_STRING_SIZE)) {
        if(productInfo->palette)
            free(productInfo->palette);
        productInfo->palette = strdup(tmpStr);
    }

    xmlGetChildStr(productNode, "category", tmpStr, XML_STRING_SIZE);
    if(productInfo->category)
        free(productInfo->category);
    productInfo->category = strdup(tmpStr);

    if(xmlGetChildStrOpt(productNode, "displayScale", tmpStr, XML_STRING_SIZE)) {
        if(productInfo->displayScale)
            free(productInfo->displayScale);
        productInfo->displayScale = strdup(tmpStr);
    }

    if(xmlGetChildStrOpt(productNode, "type", tmpStr, XML_STRING_SIZE)) {
        if(productInfo->dataType)
            free(productInfo->dataType);
        productInfo->dataType = strdup(tmpStr);
    }

    if(xmlGetChildStrOpt(productNode, "reference", tmpStr, XML_STRING_SIZE)) {
        if(productInfo->reference)
            free(productInfo->reference);
        productInfo->reference = strdup(tmpStr);
    }
    if(xmlGetChildStrOpt(productNode, "comment", tmpStr, XML_STRING_SIZE)) {
            if(productInfo->comment)
                free(productInfo->comment);
            productInfo->comment = strdup(tmpStr);
        }

    readParamDesignator(productInfo, productNode);
    if(strcmp(productInfo->paramDesignator, "none") != 0) {
        productInfo->prod_ix = paramVal;
    }

    readRange(productInfo, productNode, paramVal);
}

/**
 * read algorithm header from the XML algorithmNode and productNode
 * into algorithmInfo
 */
void readAlgorithmHeader() {
    char tmpStr[XML_STRING_SIZE];

    // copy header info from the productInfo
    copyProductInfoHeader(algorithmInfo, productInfo);

    // clear algorithm header
    if(algorithmInfo->algorithmName) {
        free(algorithmInfo->algorithmName);
        algorithmInfo->algorithmName = NULL;
    }
    if(algorithmInfo->prefix) {
        free(algorithmInfo->prefix);
        algorithmInfo->prefix = NULL;
    }
    if(algorithmInfo->suffix) {
        free(algorithmInfo->suffix);
        algorithmInfo->suffix = NULL;
    }

    // populate the fields
    if(xmlGetAttributeStrOpt(algorithmNode, "name", tmpStr, XML_STRING_SIZE)) {
        algorithmInfo->algorithmName = strdup(tmpStr);
    } else {
        algorithmInfo->algorithmName = strdup("");
    }

    readParamDesignator(algorithmInfo, algorithmNode);

    char defaultPrefix[XML_STRING_SIZE];
    if(strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        strcpy(defaultPrefix, algorithmInfo->productName);
    } else {
        strcpy(defaultPrefix, algorithmInfo->productName);
        strcat(defaultPrefix, "_");
    }

    if(xmlGetChildStrOpt(algorithmNode, "prefix", tmpStr, XML_STRING_SIZE)) {
        algorithmInfo->prefix = strdup(tmpStr);
    } else {
        algorithmInfo->prefix = strdup(defaultPrefix);
    }
    if(xmlGetChildStrOpt(algorithmNode, "suffix", tmpStr, XML_STRING_SIZE)) {
        algorithmInfo->suffix = strdup(tmpStr);
    } else {
        if(strlen(algorithmInfo->algorithmName) > 0) {
            algorithmInfo->suffix = (char*) malloc(strlen(algorithmInfo->algorithmName)+2);
            strcpy(algorithmInfo->suffix, "_");
            strcat(algorithmInfo->suffix, algorithmInfo->algorithmName);
        } else {
            algorithmInfo->suffix = strdup("");
        }
    }
}

/**
 * read algorithm data from the XML algorithmNode into algorithmInfo
 */
void readAlgorithm(int paramVal) {
    char tmpStr[XML_STRING_SIZE];

    // start with the productInfo as default values
    copyProductInfo(algorithmInfo, productInfo);
    readAlgorithmHeader();

    algorithmInfo->cat_ix = xmlGetChildDouble(algorithmNode, "cat_ix");
    xmlGetChildIntOpt(algorithmNode, "rank", &algorithmInfo->rank);

    if(xmlGetChildStrOpt(algorithmNode, "units", tmpStr, XML_STRING_SIZE)) {
        if(algorithmInfo->units)
            free(algorithmInfo->units);
        algorithmInfo->units = strdup(tmpStr);
    }

    xmlGetChildDoubleOpt(algorithmNode, "fillValue", &algorithmInfo->fillValue);

    xmlGetChildStr(algorithmNode, "description", tmpStr, XML_STRING_SIZE);
    if(algorithmInfo->description)
        free(algorithmInfo->description);
    if(strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        xmlGetChildIntOpt(algorithmNode, "prod_ix", &algorithmInfo->prod_ix);
        algorithmInfo->description = strdup(tmpStr);
    } else {
        algorithmInfo->description = (char*) malloc(strlen(tmpStr) + 64);
        algorithmInfo->prod_ix = paramVal;
        sprintf(algorithmInfo->description, tmpStr, algorithmInfo->prod_ix);
    }

    if(xmlGetChildStrOpt(algorithmNode, "reference", tmpStr, XML_STRING_SIZE)) {
        if(algorithmInfo->reference)
            free(algorithmInfo->reference);
        algorithmInfo->reference = strdup(tmpStr);
    }
    if(xmlGetChildStrOpt(algorithmNode, "comment", tmpStr, XML_STRING_SIZE)) {
            if(algorithmInfo->comment)
                free(algorithmInfo->comment);
            algorithmInfo->comment = strdup(tmpStr);
        }

    if(xmlGetChildStrOpt(algorithmNode, "type", tmpStr, XML_STRING_SIZE)) {
        if(algorithmInfo->dataType)
            free(algorithmInfo->dataType);
        algorithmInfo->dataType = strdup(tmpStr);
    }

    readRange(algorithmInfo, algorithmNode, paramVal);
}

/**
 * find first algorithm node header
 */
void findFirstAlgorithmHeader() {
    char tmpStr[XML_STRING_SIZE];

    productNode = roxml_get_chld(productsNode, NULL, 0);
    if(productNode == NULL) {
        printf("%s Line %d: could not find first product tag in XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        exit(1);
    }
    checkProductNode();
    readProductHeader();

    algorithmNode = roxml_get_chld(productNode, "algorithm", 0);
    if(algorithmNode == NULL) {
        printf("%s Line %d: could not find first algorithm tag for product=%s in XML file = %s\n",
                __FILE__, __LINE__, productInfo->productName, productXMLFileName);
        exit(1);
    }
    checkAlgorithmNode();
    readAlgorithmHeader();
}

/**
 * find next algorithm node header
 *
 * @return 1 if algorithm found, 0 if reached end of file
 */
int findNextAlgorithmHeader() {
    char tmpStr[XML_STRING_SIZE];

    algorithmNode = roxml_get_next_sibling(algorithmNode);
    if(algorithmNode == NULL) {

        // try next product
        productNode = roxml_get_next_sibling(productNode);
        if(productNode == NULL)
            return 0;
        checkProductNode();
        readProductHeader();

        algorithmNode = roxml_get_chld(productNode, "algorithm", 0);
        if(algorithmNode == NULL) {
            printf("%s Line %d: could not find first algorithm tag for product=%s in XML file = %s\n",
                    __FILE__, __LINE__, productInfo->productName, productXMLFileName);
            exit(1);
        }
    }
    checkAlgorithmNode();
    readAlgorithmHeader();
    return 1;
}

/**
 * compare algorithmInfo header to product string
 *
 * @param productFullName full name of the product to compare
 * @param sensorId sensor ID to use for wavelength comparison
 * @param paramVal pointer to an int where parameter value will be written
 * @return 1 if matches, 0 if not
 */
int compareAlgorithmHeader(const char* productFullName, int sensorId, int* paramVal) {

    // see if prefix matches
    if(strncmp(productFullName, algorithmInfo->prefix, strlen(algorithmInfo->prefix)))
        return 0;

    char tmpStr[XML_STRING_SIZE];
    if(strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        strcpy(tmpStr, algorithmInfo->prefix);
        strcat(tmpStr, algorithmInfo->suffix);
        if(strcmp(tmpStr, productFullName) == 0) {
            *paramVal = PRODUCT_DEFAULT_prod_ix;
            return 1;
        } else {
            return 0;
        }
    } else {
        int nameLen = strlen(productFullName);
        int prefixLen = strlen(algorithmInfo->prefix);
        int suffixLen = strlen(algorithmInfo->suffix);
        int paramLen = nameLen - prefixLen - suffixLen;

        // need at least one char to be a valid param
        if(paramLen < 1)
            return 0;

        // check the suffix
        if(strcmp(algorithmInfo->suffix, productFullName+nameLen-suffixLen) != 0)
            return 0;

        // get the param value
        char paramStr[XML_STRING_SIZE];
        if(paramLen>XML_STRING_SIZE-1)     // just in case
            paramLen = XML_STRING_SIZE-1;
        paramStr[0] = 0;
        strncat(paramStr, productFullName+prefixLen, paramLen);

        if(!isValidInt(paramStr))
            return 0;

        int i = atoi(paramStr);
        if(strcmp(algorithmInfo->paramDesignator, "wave") == 0) {
            if(algorithmInfo->paramWaveMin <= i && i <= algorithmInfo->paramWaveMax) {
                int32_t numBands;
                int32_t *iwave;
                int waveIndex;
                int old_verbose = want_verbose;
                want_verbose = 0;
                numBands = rdsensorinfo(sensorId, 0, "iwave", (void**)&iwave);
                numBands += rdsensorinfo(sensorId, 0, "NbandsIR", NULL);
                want_verbose = old_verbose;
                if(numBands != -1) {
                    for(waveIndex=0; waveIndex<numBands; waveIndex++) {
                        if(i == iwave[waveIndex]) {
                            *paramVal = i;
                            return 1;
                        }
                    }
                }

            }
        } else {
            *paramVal = i;
            return 1;
        }
    }
    return 0;
}

/**
 * find first product and fill in the product structure.
 *
 * @param info pre allocated product structure to fill in
 */
void getFirstProductInfo(productInfo_t* info) {
    initXmlFile();
    findFirstAlgorithmHeader();
    readProduct(-1);
    readAlgorithm(-1);
    copyProductInfo(info, algorithmInfo);
}

/**
 * find next product and fill in the product structure.
 *
 * @param info pre allocated product structure to fill in
 * @return 1 if product found, 0 if no more products
 */
int getNextProductInfo(productInfo_t* info) {
    if(findNextAlgorithmHeader()) {
        readProduct(-1);
        readAlgorithm(-1);
        copyProductInfo(info, algorithmInfo);
        return 1;
    }
    return 0;
}

/**
 * find product and fill in the product structure.
 *
 * @param productName product to find
 * @param sensorId sensor ID to use when looking up the product
 * @param info pre allocated product structure to fill (allocateProductInfo)
 * @return 1 if product found, 0 if not found
 */
int findProductInfo(const char* productName, int sensorId, productInfo_t* info) {
    int paramVal;
    const char* tmpProductName = productName;


    if(strcmp(productName, "chl_ocx") == 0) {
        switch(sensorId) {
        case OCTS:
        case SEAWIFS:
        case OCM1:
        case OCM2:
        case MOS:
        case MERIS:
        case HICO:
        case ORCA:
        case AVIRIS:
        case PRISM:
        case OLCI:
            tmpProductName = "chl_oc4";
            break;
        case HMODIST:
        case HMODISA:
        case CZCS:
        case OSMI:
        case VIIRS:
        case OCRVC:
        case GOCI:
        case OLI:
            tmpProductName = "chl_oc3";
            break;
        default:
            printf("%s Line %d: need a default chlorophyll algorithm for this sensor\n",
                    __FILE__, __LINE__);
            exit(1);
            break;

        }
    }

    initXmlFile();
    findFirstAlgorithmHeader();
    do {
        if(compareAlgorithmHeader(tmpProductName, sensorId, &paramVal)) {
            readProduct(paramVal);
            readAlgorithm(paramVal);
            copyProductInfo(info, algorithmInfo);
            return 1;
        }
    } while(findNextAlgorithmHeader());

    return 0;
}

/**
 * 
 * @param info product structure to read
 * @return full product name.  Pointer to internal memory.
 */
char* getProductNameFull(productInfo_t* info) {
    static char name[XML_STRING_SIZE];
    
    if(strcmp(info->paramDesignator, "none") == 0) {
        strcpy(name, info->prefix);
        strcat(name, info->suffix);
    } else {
        sprintf(name, "%s%d%s", info->prefix, info->prod_ix, info->suffix);
    }
    return name;
}

/**
 * print product structure.
 *
 * @param productFullName product name
 * @param info structure to print
 */
void printProductInfo(const char* productFullName, const productInfo_t* info) {
    printf("fullProductName=%s\n", productFullName);
    printf("description=%s\n", info->description);
    printf("units=%s\n", info->units);
    printf("palette=%s\n", info->palette);
    printf("paramDesignator=%s\n", info->paramDesignator);
    printf("paramWaveMin=%d\n", info->paramWaveMin);
    printf("paramWaveMax=%d\n", info->paramWaveMax);
    printf("standardName=%s\n", info->standardName);
    printf("category=%s\n", info->category);
    printf("dataType=%s\n", info->dataType);
    printf("prefix=%s\n", info->prefix);
    printf("suffix=%s\n", info->suffix);
    printf("algorithmName=%s\n", info->algorithmName);
    printf("productName=%s\n", info->productName);
    printf("cat_ix=%d\n", info->cat_ix);
    printf("prod_ix=%d\n", info->prod_ix);
    printf("rank=%d\n", info->rank);
    printf("fillValue=%g\n", info->fillValue);
    printf("validMin=%g\n", info->validMin);
    printf("validMax=%g\n", info->validMax);
    printf("displayScale=%s\n", info->displayScale);
    printf("displayMin=%g\n", info->displayMin);
    printf("displayMax=%g\n", info->displayMax);
    printf("scaleFactor=%g\n", info->scaleFactor);
    printf("reference=%s\n", info->reference);
    printf("comment=%s\n", info->comment);
}

