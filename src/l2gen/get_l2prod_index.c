/* -------------------------------------------------------------------------- */
/* Module get_l2prod_index - maintains a list of all possible L2 products and */
/* returns attributes on request.                                             */
/*                                                                            */
/* B. A. Franz, SAIC GSC, NASA/SIMBIOS Project, April 1999                    */
/* -------------------------------------------------------------------------- */

#include <roxml.h>
#include <xmlUtils.h>

#include <stdio.h>
#include <strings.h>
#include "l12_proto.h"
#include "l2prod.h"



#define ARRAY_CHUNK_SIZE 20

static l2prodstr **l2prod_array = NULL; // array of pointers to all the products
static int l2prod_num = 0; // number of products in the array
static int l2prod_storage = 0; // size of the array holding the pointers

/* init a product structure */
void initProduct(l2prodstr* l2prod) {
    l2prod->param_type = PARAM_TYPE_NONE;
    l2prod->name_prefix[0] = '\0';
    l2prod->name_suffix[0] = '\0';
    l2prod->cat_ix = -1;
    l2prod->prod_ix = -1;
    l2prod->datatype = DFNT_FLOAT32;
    l2prod->slope = 1.0;
    l2prod->offset = 0.0;
    l2prod->min = 0;
    l2prod->max = 0;
    l2prod->rank = 2;
    strcpy(l2prod->title_format, "no title format yet (%d)");
    strcpy(l2prod->title, "no title yet");
    strcpy(l2prod->units, "undefined units");
    l2prod->badData = BAD_FLT;
    l2prod->product_id[0] = '\0';
    l2prod->algorithm_id[0] = '\0';
    l2prod->standard_name[0] = '\0';
}

/* create a new product structure, init it, add it to the global array
 * and return the pointer.
 */
l2prodstr* createNewProduct() {

    // allocate the structure
    l2prodstr* l2prod = (l2prodstr*) malloc(sizeof (l2prodstr));

    initProduct(l2prod);
    
    // add the product structure to the global array
    l2prod_num++;
    if (l2prod_num > l2prod_storage) {
        l2prod_storage += ARRAY_CHUNK_SIZE;
        l2prod_array = (l2prodstr**) realloc(l2prod_array, l2prod_storage * sizeof (l2prodstr*));
    }
    l2prod_array[l2prod_num - 1] = l2prod;

    return l2prod;
}


/* -------------------------------------------------------------------------- */
int convertDataType(char* str) {
    if(strcmp(str, "byte") == 0)
        return DFNT_INT8;
    if(strcmp(str, "ubyte") == 0)
        return DFNT_UINT8;
    if(strcmp(str, "short") == 0)
        return DFNT_INT16;
    if(strcmp(str, "ushort") == 0)
        return DFNT_UINT16;
    if(strcmp(str, "int") == 0)
        return DFNT_INT32;
    if(strcmp(str, "uint") == 0)
        return DFNT_UINT32;
    if(strcmp(str, "float") == 0)
        return DFNT_FLOAT32;
    if(strcmp(str, "double") == 0)
        return DFNT_FLOAT64;
    return DFNT_FLOAT32;
}

/* -------------------------------------------------------------------------- */
int convertParamType(char* str) {
    if(strcmp(str, "none") == 0)
        return PARAM_TYPE_NONE;
    if(strcmp(str, "uv") == 0)
        return PARAM_TYPE_ALL_WAVE;
    if(strcmp(str, "visible") == 0)
        return PARAM_TYPE_VIS_WAVE;
    if(strcmp(str, "nir") == 0)
        return PARAM_TYPE_ALL_WAVE;
    if(strcmp(str, "swir") == 0)
        return PARAM_TYPE_ALL_WAVE;
    if(strcmp(str, "emissive") == 0)
        return PARAM_TYPE_ALL_WAVE;
    if(strcmp(str, "band") == 0)
        return PARAM_TYPE_BAND;
    if(strcmp(str, "int") == 0)
        return PARAM_TYPE_INT;
    return PARAM_TYPE_NONE;
}

/* -------------------------------------------------------------------------- */
void init_l2prod() {
    int i;
    l2prodstr *l2prod;
    node_t* rootNode;
    node_t* productsNode;
    node_t* productNode;
    node_t* algorithmNode;
    node_t* rangeNode;
    node_t* tmpNode;
    char productName[UNITLEN];
    char standardName[TITLELEN];
    char productUnits[UNITLEN];
    char productDataType[UNITLEN];
    char productParamDesignator[UNITLEN];
    double productMin, productMax;
    double productScale, productOffset;
    char algorithmName[UNITLEN];
    char tmpStr[TITLELEN];
    char fileName[FILENAME_MAX];

    char *dataRoot;
    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }
    strcpy(fileName, dataRoot);
    strcat(fileName,"/common/product.xml");

    rootNode = roxml_load_doc(fileName);
    if(rootNode == NULL) {
       printf("%s Line %d: could not open product XML file = %s\n",
               __FILE__, __LINE__, fileName);
       exit(1);
    }

    productsNode = roxml_get_chld(rootNode, "products", 0);
    productNode = roxml_get_chld(productsNode, NULL, 0);
    while(productNode) {
        xmlGetAttributeStr(productNode, "name", productName, UNITLEN);
        standardName[0] = 0;
        xmlGetChildStrOpt(productNode, "standardName", standardName, TITLELEN);
        xmlGetChildStr(productNode, "units", productUnits, TITLELEN);
        strcpy(productDataType, "float");
        xmlGetChildStrOpt(productNode, "type", productDataType, UNITLEN);
        strcpy(productParamDesignator, "none");
        xmlGetChildStrOpt(productNode, "paramDesignator", productParamDesignator, UNITLEN);
        productMin = 0;
        productMax = 0;
        productScale = 1;
        productOffset = 0;

        rangeNode = roxml_get_chld(productNode, "range", 0);
        if(rangeNode) {
            xmlGetChildDoubleOpt(rangeNode, "validMin", &productMin);
            xmlGetChildDoubleOpt(rangeNode, "validMax", &productMax);
            xmlGetChildDoubleOpt(rangeNode, "scaleFactor", &productScale);
            xmlGetChildDoubleOpt(rangeNode, "addOffset", &productOffset);
        }

        algorithmNode = roxml_get_chld(productNode, "algorithm", 0);
        while(algorithmNode != NULL) {
            roxml_get_name(algorithmNode, tmpStr, TITLELEN);
            if(strcmp(tmpStr, "algorithm") != 0) {
                printf("%s Line %d: \"algorithm\" expected, found \"%s\"\n",
                        __FILE__, __LINE__, tmpStr);
                exit(1);
            }

            l2prod = createNewProduct();
            algorithmName[0] = 0;
            xmlGetAttributeStrOpt(algorithmNode, "name", algorithmName, UNITLEN);
            strcpy(l2prod->product_id, productName);
            strcpy(l2prod->algorithm_id, algorithmName);
            strcpy(l2prod->standard_name, standardName);
            l2prod->cat_ix = xmlGetChildDouble(algorithmNode, "cat_ix");
            xmlGetChildIntOpt(algorithmNode, "prod_ix", &l2prod->prod_ix);
            xmlGetChildIntOpt(algorithmNode, "rank", &l2prod->rank);
            strcpy(l2prod->units, productUnits);
            xmlGetChildStrOpt(algorithmNode, "units", l2prod->units, UNITLEN);
            xmlGetChildStr(algorithmNode, "description", l2prod->title_format, TITLELEN);
            l2prod->badData = BAD_FLT;
            xmlGetChildFloatOpt(algorithmNode, "fillValue", &l2prod->badData);
            strcpy(tmpStr, productDataType);
            xmlGetChildStrOpt(algorithmNode, "type", tmpStr, TITLELEN);
            l2prod->datatype = convertDataType(tmpStr);
            strcpy(tmpStr, productParamDesignator);
            xmlGetChildStrOpt(algorithmNode, "paramDesignator", tmpStr, TITLELEN);
            l2prod->param_type = convertParamType(tmpStr);

            // set values that depend on param_type
            int prefixExists = xmlGetChildStrOpt(algorithmNode, "prefix", l2prod->name_prefix, UNITLEN);
            int suffixExists = xmlGetChildStrOpt(algorithmNode, "suffix", l2prod->name_suffix, UNITLEN);
            if(l2prod->param_type == PARAM_TYPE_NONE) {
                if(!prefixExists) {
                    strcpy(l2prod->name_prefix, productName);
                }
                if(!suffixExists && algorithmName[0] != 0) {
                    strcpy(l2prod->name_suffix, "_");
                    strcat(l2prod->name_suffix, algorithmName);
                }
                strcpy(l2prod->title, l2prod->title_format);
            } else {
                if(!prefixExists) {
                    strcpy(l2prod->name_prefix, productName);
                    strcat(l2prod->name_prefix, "_");
                }
                if(!suffixExists && algorithmName[0] != 0) {
                    strcpy(l2prod->name_suffix, "_");
                    strcat(l2prod->name_suffix, algorithmName);
                }
            }

            l2prod->min = productMin;
            l2prod->max = productMax;
            l2prod->slope = productScale;
            l2prod->offset = productOffset;
            rangeNode = roxml_get_chld(algorithmNode, "range", 0);
            if(rangeNode) {
                xmlGetChildDoubleOpt(rangeNode, "validMin", &l2prod->min);
                xmlGetChildDoubleOpt(rangeNode, "validMax", &l2prod->max);
                xmlGetChildFloatOpt(rangeNode, "scaleFactor", &l2prod->slope);
                xmlGetChildFloatOpt(rangeNode, "addOffset", &l2prod->offset);
            }

            algorithmNode = roxml_get_next_sibling(algorithmNode);
        }

        // go to the next product
        productNode = roxml_get_next_sibling(productNode);
    }

    roxml_close(rootNode);
}


/* read the product structure and copy the product name into name */
void getProductName(char* name, l2prodstr* product) {
    if (product->product_id[0] != '\0') {
        strcpy(name, product->product_id);
    } else {
        strcpy(name, product->name_prefix);

        // remove trailing underscore
        if (name[strlen(name) - 1] == '_') {
            name[strlen(name) - 1] = '\0';
        }
    }
}

/* read the product structure and copy the prefix into name */
void getPrefix(char* name, l2prodstr* product) {
    char productName[UNITLEN];
    getProductName(productName, product);
    strcat(productName, "_");
    if(strcmp(productName, product->name_prefix) == 0)
        name[0] = 0;
    else
        strcpy(name, product->name_prefix);
    return;
}

/* read the product structure and copy the algorithm name into name */
void getAlgorithmName(char* name, l2prodstr* product) {
    int i;
    int length;

    if (product->algorithm_id[0] != '\0') {
        strcpy(name, product->algorithm_id);
    } else {
        strcpy(name, product->name_suffix);

        // remove leading underscore
        if (name[0] == '_') {
            length = strlen(name);
            for (i = 1; i <= length; i++) {
                name[i - 1] = name[i];
            }
        }
    }
}

/* read the product structure and copy the suffix into name */
void getSuffix(char* name, l2prodstr* product) {
    char algorithmName[UNITLEN];
    algorithmName[0] = '_';
    getAlgorithmName(algorithmName+1, product);
    if(strcmp(algorithmName, product->name_suffix) == 0)
        name[0] = 0;
    else
        strcpy(name, product->name_suffix);
    return;
}

/* read the product structure and copy the product param type into name */
void getParamType(char* name, l2prodstr* product) {
    switch (product->param_type) {
        case PARAM_TYPE_NONE:
            strcpy(name, "None");
            break;
        case PARAM_TYPE_VIS_WAVE:
            strcpy(name, "Visible");
            break;
        case PARAM_TYPE_ALL_WAVE:
            strcpy(name, "All");
            break;
        case PARAM_TYPE_BAND:
            strcpy(name, "Bands");
            break;
        case PARAM_TYPE_INT:
            strcpy(name, "Integer");
            break;
        default:
            strcpy(name, "Unknown");
            break;
    }
}

/* read the product structure and copy the product data type into name */
void getDataType(char* name, l2prodstr* product) {
    switch (product->datatype) {
        case DFNT_INT8:
            strcpy(name, "INT8");
            break;
        case DFNT_UINT8:
            strcpy(name, "UINT8");
            break;
        case DFNT_INT16:
            strcpy(name, "INT16");
            break;
        case DFNT_INT32:
            strcpy(name, "INT32");
            break;
        case DFNT_FLOAT32:
            strcpy(name, "FLOAT32");
            break;
        case DFNT_FLOAT64:
            strcpy(name, "FLOAT64");
            break;
        default:
            strcpy(name, "Unknown");
            break;
    }
}

void setTitleString(l2prodstr* product) {
    switch (product->param_type) {
        case PARAM_TYPE_NONE:
            break;
        case PARAM_TYPE_VIS_WAVE:
        case PARAM_TYPE_ALL_WAVE:
        case PARAM_TYPE_BAND:
        case PARAM_TYPE_INT:
            sprintf(product->title, product->title_format, 0);
            break;
    }
}

void getFullProductName(char* name, l2prodstr* prod) {
    if(prod->param_type == PARAM_TYPE_NONE) {
        strcpy(name, prod->name_prefix);
        strcat(name, prod->name_suffix);
    } else {
        strcpy(name, prod->name_prefix);
        strcat(name, "nnn");
        strcat(name, prod->name_suffix);
    }
}

/* dump the product list to a file */
void dumpProductStructure(l2prodstr **list, char* filename) {
    l2prodstr* product;
    int productIndex;
    char fullProductName[UNITLEN];
    char productName[UNITLEN];
    char algorithmName[UNITLEN];
    char paramType[UNITLEN];
    char* description;
    FILE* fout;

    fout = fopen(filename, "w");

    for (productIndex = 0; productIndex < l2prod_num; productIndex++) {
        product = list[productIndex];
        getProductName(productName, product);
        getAlgorithmName(algorithmName, product);
        getParamType(paramType, product);
        getFullProductName(fullProductName, product);
        if(product->param_type == PARAM_TYPE_NONE) {
            description = product->title;
        } else {
            description =  product->title_format;
        }

        fprintf(fout, "%s, %s, %s, %s, %d, %d, %d, %g, %g, %g, %g, %g, %d, %s, %s, %s\n",
        fullProductName,
        paramType,
        productName,
        algorithmName,
        product->cat_ix,
        product->prod_ix,
        product->datatype,
        product->slope,
        product->offset,
        product->min,
        product->max,
        product->badData,
        product->rank,
        product->units,
        description,
        product->standard_name);
    }

    fclose(fout);
}

char* xmlEncapsulateText(char* str) {
    static int strLength = 0;
    static char* buffer = NULL;

    if (strpbrk(str, "<>&")) {
        // make buffer big enough
        int length = strlen(str) + 15;
        if(strLength < length) {
            if(buffer)
                free(buffer);
            strLength = length;
            buffer = (char*) malloc(strLength);
        }

        sprintf(buffer, "<![CDATA[%s]]>", str);
        return buffer;
    }
    return str;
}

void XML_file_header(FILE* fout) {
    fprintf(fout, "<productList>\n");
}

void XML_file_footer(FILE* fout) {
    fprintf(fout, "</productList>\n");
}

void XML_product_header(FILE* fout, l2prodstr* product) {
    char name[UNITLEN];

    getProductName(name, product);
    fprintf(fout, "    <product name=\"%s\">\n", name);
}

void XML_product_footer(FILE* fout) {
    fprintf(fout, "    </product>\n");
}

void XML_algorithm(FILE* fout, l2prodstr* product) {
    char name[UNITLEN];
    char param[UNITLEN];
    char data[UNITLEN];

    getAlgorithmName(name, product);
    getParamType(param, product);
    getDataType(data, product);
    setTitleString(product);

    if (strlen(name) > 0) {
        fprintf(fout, "        <algorithm name=\"%s\">\n", name);
    } else {
        fprintf(fout, "        <algorithm>\n");
    }
    clo_writeXmlTag(fout, 3, "prefix", product->name_prefix);
    clo_writeXmlTag(fout, 3, "suffix", product->name_suffix);
    clo_writeXmlTag(fout, 3, "parameterType", param);
    clo_writeXmlTag(fout, 3, "dataType", data);
    clo_writeXmlTag(fout, 3, "units", product->units);
    clo_writeXmlTag(fout, 3, "description", product->title);
    fprintf(fout, "        </algorithm>\n");
}

void write_product_XML_file(char* filename) {
    l2prodstr **tmpList = NULL; // array of pointers to all the products
    l2prodstr **sortedList = NULL; // array of pointers to all the products
    int sortedIndex = 0; // number of products in the array
    l2prodstr* product;
    l2prodstr* product2;
    static l2prodstr *chlProd;
    char productName[UNITLEN];
    char productName2[UNITLEN];
    int productIndex;
    int productIndex2;
    FILE* fout;
    int i;

    if(filename == NULL || filename[0] == 0) {
        printf("%s Line %d: NULL filename for XML file.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if(chlProd==NULL) {
        // allocate the structure
        chlProd = (l2prodstr*) malloc(sizeof (l2prodstr));
        initProduct(chlProd);

        chlProd->cat_ix = CAT_chl_oc2;
        strncpy(chlProd->name_prefix, "chlor_a", UNITLEN);
        strncpy(chlProd->product_id, "chlor_a", UNITLEN);
        strncpy(chlProd->algorithm_id, "chlor_a", UNITLEN);
        strncpy(chlProd->units, "mg m^-3", UNITLEN);
        strncpy(chlProd->title, "Chlorophyll Concentration, Default Sensor Algorithm", TITLELEN);
        strncpy(chlProd->standard_name, "chlorophyll_concentration_in_sea_water", TITLELEN);
    }

    // allocate lists
    tmpList = (l2prodstr**) malloc((l2prod_num+1) * sizeof (l2prodstr*));
    sortedList = (l2prodstr**) malloc((l2prod_num+1) * sizeof (l2prodstr*));

    // copy list
    for (productIndex = 0; productIndex < l2prod_num; productIndex++) {
        tmpList[productIndex] = l2prod_array[productIndex];
    }

    // add default chl
    tmpList[productIndex] = chlProd;

    // sort list by product name
    for (productIndex = 0; productIndex < l2prod_num+1; productIndex++) {
        product = tmpList[productIndex];
        if (product) {
            getProductName(productName, product);
            sortedList[sortedIndex++] = product; // add to sorted
            tmpList[productIndex] = NULL; // remove from tmp

            // search for matching product entries
            for (productIndex2 = productIndex + 1; productIndex2 < l2prod_num+1; productIndex2++) {
                product2 = tmpList[productIndex2];
                if (product2) {
                    getProductName(productName2, product2);
                    if (strcasecmp(productName, productName2) == 0) {
                        sortedList[sortedIndex++] = product2; // add to sorted
                        tmpList[productIndex2] = NULL; // remove from tmp
                    } // name matches
                } // product2 exists
            } // for product2
        } // product exists
    } // for product

    // open XML file
    fout = fopen(filename, "w");

    // loop sorted list to print XML
    XML_file_header(fout);
    getProductName(productName2, sortedList[0]);
    XML_product_header(fout, sortedList[0]);

    for (productIndex = 0; productIndex < l2prod_num+1; productIndex++) {
        product = sortedList[productIndex];
        getProductName(productName, product);

        // if changed product
        if (strcasecmp(productName, productName2) != 0) {
            strcpy(productName2, productName);
            XML_product_footer(fout);
            XML_product_header(fout, product);
        }
        XML_algorithm(fout, product);
    }

    XML_product_footer(fout);
    XML_file_footer(fout);

    // close XML file
    fclose(fout);

    // free lists
    free(tmpList);
    free(sortedList);
}

node_t* xmlMakeProduct(node_t* rootNode, char* productName, l2prodstr* product) {
    node_t* productNode = roxml_add_node(rootNode, 0, ROXML_ELM_NODE, "product", NULL);
    roxml_add_node(productNode, 0, ROXML_ATTR_NODE, "name", productName);
    if(product->standard_name[0] != 0)
        roxml_add_node(productNode, 0, ROXML_ELM_NODE, "standardName", product->standard_name);
    roxml_add_node(productNode, 0, ROXML_ELM_NODE, "palette", "defaultPalette");
    roxml_add_node(productNode, 0, ROXML_ELM_NODE, "category", "defaultCategory");
    return productNode;
}

node_t* xmlMakeAlgorithm(node_t* productNode, l2prodstr* product) {
    char nameStr[UNITLEN];
    char tmpStr[2048];

    node_t* algorithmNode = roxml_add_node(productNode, 0, ROXML_ELM_NODE, "algorithm", NULL);
    getAlgorithmName(nameStr, product);
    if(nameStr[0] != 0)
        roxml_add_node(algorithmNode, 0, ROXML_ATTR_NODE, "name", nameStr);

    sprintf(tmpStr, "%d", product->cat_ix);
    roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "cat_ix", tmpStr);
    if(product->prod_ix != -1) {
        sprintf(tmpStr, "%d", product->prod_ix);
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "prod_ix", tmpStr);
    }
    if(product->rank != 2) {
        sprintf(tmpStr, "%d", product->rank);
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "rank", tmpStr);
    }
    getPrefix(nameStr, product);
    if(nameStr[0] != 0) {
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "prefix", nameStr);
    }
    getSuffix(nameStr, product);
    if(nameStr[0] != 0) {
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "suffix", nameStr);
    }
    roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "units", xmlEncapsulateText(product->units));
    sprintf(tmpStr, "%g", product->badData);
    roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "fillValue", tmpStr);
    if(product->param_type == PARAM_TYPE_NONE) {
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "description", xmlEncapsulateText(product->title));
    } else {
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "description", xmlEncapsulateText(product->title_format));
    }

    if(product->datatype != DFNT_FLOAT32) {
        getDataType(nameStr, product);
        roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "type", nameStr);
    }

    if(product->param_type != PARAM_TYPE_NONE) {
        switch(product->param_type) {
        case PARAM_TYPE_VIS_WAVE:
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "visible");
           break;
        case PARAM_TYPE_ALL_WAVE:
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "uv");
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "visible");
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "nir");
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "swir");
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "emissive");
            break;
        case PARAM_TYPE_BAND:
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "band");
            break;
        case PARAM_TYPE_INT:
            roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "paramDesignator", "int");
            break;
        default:
            printf("Error - paramType not recognized\n");
            exit(1);
        }
    }

    node_t* rangeNode = roxml_add_node(algorithmNode, 0, ROXML_ELM_NODE, "range", NULL);
    sprintf(tmpStr, "%g", product->min);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "validMin", tmpStr);
    sprintf(tmpStr, "%g", product->max);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "validMax", tmpStr);
    sprintf(tmpStr, "%g", product->min);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "displayMin", tmpStr);
    sprintf(tmpStr, "%g", product->max);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "displayMax", tmpStr);
    sprintf(tmpStr, "%g", product->slope);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "scaleFactor", tmpStr);
    sprintf(tmpStr, "%g", product->offset);
    roxml_add_node(rangeNode, 0, ROXML_ELM_NODE, "addOffset", tmpStr);

    return algorithmNode;
}

void write_product_XML_file2(char* filename) {
    l2prodstr **tmpList = NULL; // array of pointers to all the products
    l2prodstr **sortedList = NULL; // array of pointers to all the products
    int sortedIndex = 0; // number of products in the array
    l2prodstr* product;
    l2prodstr* product2;
    static l2prodstr *chlProd;
    char productName[UNITLEN];
    char productName2[UNITLEN];
    int productIndex;
    int productIndex2;
    FILE* fout;
    int i;

    if(filename == NULL || filename[0] == 0) {
        printf("%s Line %d: NULL filename for XML file.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if(chlProd==NULL) {
        // allocate the structure
        chlProd = (l2prodstr*) malloc(sizeof (l2prodstr));
        initProduct(chlProd);

        chlProd->cat_ix = CAT_chl_oc2;
        strncpy(chlProd->name_prefix, "chlor_a", UNITLEN);
        strncpy(chlProd->product_id, "chlor_a", UNITLEN);
        strncpy(chlProd->algorithm_id, "chlor_a", UNITLEN);
        strncpy(chlProd->units, "mg m^-3", UNITLEN);
        strncpy(chlProd->title, "Chlorophyll Concentration, Default Sensor Algorithm", TITLELEN);
        strncpy(chlProd->standard_name, "chlorophyll_concentration_in_sea_water", TITLELEN);
    }

    // allocate lists
    tmpList = (l2prodstr**) malloc((l2prod_num+1) * sizeof (l2prodstr*));
    sortedList = (l2prodstr**) malloc((l2prod_num+1) * sizeof (l2prodstr*));

    // copy list
    for (productIndex = 0; productIndex < l2prod_num; productIndex++) {
        tmpList[productIndex] = l2prod_array[productIndex];
    }

    // add default chl
    tmpList[productIndex] = chlProd;

    // sort list by product name
    for (productIndex = 0; productIndex < l2prod_num+1; productIndex++) {
        product = tmpList[productIndex];
        if (product) {
            getProductName(productName, product);
            sortedList[sortedIndex++] = product; // add to sorted
            tmpList[productIndex] = NULL; // remove from tmp

            // search for matching product entries
            for (productIndex2 = productIndex + 1; productIndex2 < l2prod_num+1; productIndex2++) {
                product2 = tmpList[productIndex2];
                if (product2) {
                    getProductName(productName2, product2);
                    if (strcasecmp(productName, productName2) == 0) {
                        sortedList[sortedIndex++] = product2; // add to sorted
                        tmpList[productIndex2] = NULL; // remove from tmp
                    } // name matches
                } // product2 exists
            } // for product2
        } // product exists
    } // for product

    // make XML tree
    node_t* rootNode = roxml_add_node(NULL, 0, ROXML_PI_NODE, "xml", NULL);
    roxml_add_node(rootNode, 0, ROXML_ATTR_NODE, "version", "1.0");

    rootNode = roxml_add_node(NULL, 0, ROXML_ELM_NODE, "products", NULL);
    roxml_add_node(rootNode, 0, ROXML_ATTR_NODE, "xmlns", "http://oceancolor.gsfc.nasa.gov");
    roxml_add_node(rootNode, 0, ROXML_ATTR_NODE, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    roxml_add_node(rootNode, 0, ROXML_ATTR_NODE, "xsi:schemaLocation", "http://oceancolor.gsfc.nasa.gov product-1.0.xsd");

    node_t * productNode = NULL;

    for (productIndex = 0; productIndex < l2prod_num+1; productIndex++) {
        product = sortedList[productIndex];
        getProductName(productName, product);

        // if changed product
        if (strcasecmp(productName, productName2) != 0) {
             strcpy(productName2, productName);
             productNode = xmlMakeProduct(rootNode, productName, product);
        }
        xmlMakeAlgorithm(productNode, product);

    } // products

    roxml_commit_changes(rootNode, filename, NULL, 1);
    roxml_release(RELEASE_ALL);

    // free lists
    free(tmpList);
    free(sortedList);
}


/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
l2prodstr *get_l2prod_index(char *prod_name, /* Input SDS name              */
        int32 sensorID, /* Sensor ID                   */
        int32 numBands, /* Number of Wavelengths       */
        int32 numPixels, /* Number of Pixels per Scan   */
        int32 numScans, /* Number of Scans in File     */
        int32_t wave[]) /* Wavelength Array (numBands) */ {
    int i;
    int prodIndex;
    char tmp_pname[32];
    l2prodstr *p = NULL;
    int status;

    static int first = 1;
    static int lastNumPixels = 0;
    static int default_iprod_chl = -1;
    static l2prodstr* default_prod_chl = NULL;
    static productInfo_t* info;


    if (first) {
        first = 0;
        lastNumPixels = numPixels;
        info = allocateProductInfo();
        init_l2prod();

        // loop through the products and set some default info
        for (i = 0; i < l2prod_num; i++) {
            l2prodstr *l2prod = l2prod_array[i];
            if (l2prod->rank == 2) {
                l2prod->dim[0] = numScans;
                l2prod->dim[1] = numPixels;
                l2prod->dim[2] = 1;
                strncpy(l2prod->dimname[0], "Number of Scan Lines", DIMNAMELEN);
                strncpy(l2prod->dimname[1], "Pixels per Scan Line", DIMNAMELEN);
                strncpy(l2prod->dimname[2], "", DIMNAMELEN);
            } else if (l2prod->rank == 1) {
                l2prod->dim[0] = numScans;
                l2prod->dim[1] = 1;
                l2prod->dim[2] = 1;
                strncpy(l2prod->dimname[0], "Number of Scan Lines", DIMNAMELEN);
                strncpy(l2prod->dimname[1], "", DIMNAMELEN);
                strncpy(l2prod->dimname[2], "", DIMNAMELEN);
            }
        }
// the following needs to be consistent with libl2/productInfo.c:findProductInfo
        switch (sensorID) {
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
                default_iprod_chl = CAT_chl_oc4;
                break;
            case HMODIST:
            case HMODISA:
            case CZCS:
            case OSMI:
            case VIIRS:
            case OCRVC:
            case GOCI:
            case OLI:
                default_iprod_chl = CAT_chl_oc3;
                break;
            case AVHRR:
                break;
            default:
                printf("%s Line %d: need a default chlorophyll algorithm for this sensor\n",
                        __FILE__, __LINE__);
                exit(1);
                break;
        } // switch

        for (prodIndex = 0; prodIndex < l2prod_num; prodIndex++) {
            p = l2prod_array[prodIndex];
            if (p->cat_ix == default_iprod_chl) {
                default_prod_chl = p;
                break;
            }
        }
        if (default_prod_chl == NULL && sensorID != AVHRR) {
            printf("%s Line %d: could not find the default chlorophyll algorithm for this sensor\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        // fix the title of the CDOM-corrected CHL model
        for (prodIndex = 0; prodIndex < l2prod_num; prodIndex++) {
            p = l2prod_array[prodIndex];
            if (p->cat_ix == CAT_chl_cdomcorr_morel) {
                sprintf(p->title, "%s, CDOM-corrected via Morel",
                        default_prod_chl->title);
                break;
            }
        }

    } // if first

    if (numPixels != lastNumPixels){
        lastNumPixels = numPixels;
        for (i = 0; i < l2prod_num; i++) {
            l2prodstr *l2prod = l2prod_array[i];
            if (l2prod->rank == 2) {
                l2prod->dim[1] = numPixels;
            }
        }
    }

    if (strcmp(prod_name, "chl_ocx") == 0) { /* Sensor default chlorophyll-a */
        return default_prod_chl;
    }

    // loop through products until we find a match
    for (prodIndex = 0; prodIndex < l2prod_num; prodIndex++) {
        p = l2prod_array[prodIndex];
        switch (p->param_type) {
            case PARAM_TYPE_NONE:
                sprintf(tmp_pname, "%s%s", p->name_prefix, p->name_suffix);
                if (strcmp(prod_name, tmp_pname) == 0)
                    return p;
                break;

            case PARAM_TYPE_VIS_WAVE:
            case PARAM_TYPE_ALL_WAVE:
                if(strncmp(prod_name, p->name_prefix, strlen(p->name_prefix)) == 0) {
                    for (i = 0; i < numBands; i++) {
                        sprintf(tmp_pname, "%s%d%s", p->name_prefix, wave[i], p->name_suffix);
                        if (strcmp(prod_name, tmp_pname) == 0) {
                            p->prod_ix = i;
                            sprintf(p->title, p->title_format, wave[i]);

                            // fix the values that change due to range
                            status = findProductInfo(prod_name, sensorID, info);
                            if(status) {
                                p->slope = info->scaleFactor;
                                p->offset = info->addOffset;
                                p->max = info->validMax;
                                p->min = info->validMin;
                            }
                            return p;
                        }
                    } // for bands
                } // if prefix matches
                break;

            case PARAM_TYPE_BAND:
                for (i = 0; i < numBands; i++) {
                    sprintf(tmp_pname, "%s%d%s", p->name_prefix, i+1, p->name_suffix);
                    if (strcmp(prod_name, tmp_pname) == 0) {
                        p->prod_ix = i;
                        sprintf(p->title, p->title_format, i);

                        // fix the values that change due to range
                        status = findProductInfo(prod_name, sensorID, info);
                        if(status) {
                            p->slope = info->scaleFactor;
                            p->offset = info->addOffset;
                            p->max = info->validMax;
                            p->min = info->validMin;
                        }
                        return p;
                    }
                }
                break;

            case PARAM_TYPE_INT:
                if (p->name_suffix[0] == 0) {
                    if (strncmp(prod_name, p->name_prefix, strlen(p->name_prefix)) == 0) {
                        p->prod_ix = atoi(&prod_name[strlen(p->name_prefix)]);
                        sprintf(p->title, p->title_format, p->prod_ix);
                        // fix the values that change due to range
                        status = findProductInfo(prod_name, sensorID, info);
                        if(status) {
                            p->slope = info->scaleFactor;
                            p->offset = info->addOffset;
                            p->max = info->validMax;
                            p->min = info->validMin;
                        }
                       return p;
                    }
                } else {
                    if ((strncmp(prod_name, p->name_prefix, strlen(p->name_prefix)) == 0) &&
                            (strcmp(prod_name + strlen(prod_name) - strlen(p->name_suffix), p->name_suffix) == 0)) {
                        p->prod_ix = atoi(&prod_name[strlen(p->name_prefix)]);
                        sprintf(p->title, p->title_format, p->prod_ix);
                        // fix the values that change due to range
                        status = findProductInfo(prod_name, sensorID, info);
                        if(status) {
                            p->slope = info->scaleFactor;
                            p->offset = info->addOffset;
                            p->max = info->validMax;
                            p->min = info->validMin;
                        }
                        return p;
                    }
                } // prefix and suffix
                break;

        } // switch

    } // for product list

    fprintf(stderr, "\n%s - Invalid product name : %s\n", __FILE__, prod_name);
    exit(1);
}

