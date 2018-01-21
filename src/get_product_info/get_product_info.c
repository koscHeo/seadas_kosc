#include <stdio.h>
#include <string.h>

#include <genutils.h>
#include <clo.h>
#include <sensorInfo.h>
#include <productInfo.h>

#define VERSION "2.0"
#define PROG_NAME "get_product_info"

void printProductList() {
    char tmpStr[1024];
    productInfo_t* productInfo = allocateProductInfo();

    printf("chlor_a\n");
    getFirstProductInfo(productInfo);
    do {
        if(strcmp(productInfo->paramDesignator, "none") == 0) {
            printf("%s%s\n", productInfo->prefix, productInfo->suffix);
        } else {
            printf("%snnn%s\n", productInfo->prefix, productInfo->suffix);
        }
    } while(getNextProductInfo(productInfo));
}

void printProductListSensor(int sensor) {
    char tmpStr[1024];

    // get sensor wavelengths
    int32_t* iwave;

    // stop verbose output
    int old_verbose = want_verbose;
    want_verbose = 0;
    int numWavelengths = rdsensorinfo(sensor, 0, "iwave", (void**)&iwave);
    want_verbose = old_verbose;

    if(numWavelengths == -1) {
        printf("-E- Could not lookup sensor %d wavelengths\n", sensor);
        exit(1);
    }

    printf("chlor_a\n");
    productInfo_t* productInfo = allocateProductInfo();
    getFirstProductInfo(productInfo);
    do {
        if(strcmp(productInfo->paramDesignator, "none") == 0) {
            printf("%s%s\n", productInfo->prefix, productInfo->suffix);
        } else if(strcmp(productInfo->paramDesignator, "wave") == 0) {
            int i;
            for(i=0; i<numWavelengths; i++) {
                printf("%s%d%s\n", productInfo->prefix, iwave[i], productInfo->suffix);
            }
        } else {
            printf("%snnn%s\n", productInfo->prefix, productInfo->suffix);
        }
    } while(getNextProductInfo(productInfo));
}


int main(int argc, char *argv[])
{
    // setup clo
    clo_setVersion2(PROG_NAME, VERSION);
    char helpStr[2048];
    strcpy(helpStr, "Program to list products or detailed information about a single product\n");
    strcat(helpStr, "Usage: ");
    strcat(helpStr, PROG_NAME);
    strcat(helpStr, " [option=val] <productName>\nOptions:");
    clo_setHelpStr(helpStr);
    clo_optionList_t* list = clo_createList();

    clo_addOption(list, "-l", CLO_TYPE_BOOL, "false", "list all of the products");
    clo_addOption(list, "-r", CLO_TYPE_BOOL, "false", "list all of the products recursing through\n        the wavelengths of the sensor specified");
    clo_addOption(list, "sensor", CLO_TYPE_STRING, "MODISA", "sensor name (or ID) to use for wavelength\n        expansion or product lookup");
    clo_addOption(list, "<product>", CLO_TYPE_HELP, NULL, "<productName> product to print detailed information about");

    int initialNumOptions = clo_getNumOptions(list);
    clo_readArgs(list, argc, argv);
    int numExtraOptions = clo_getNumOptions(list) - initialNumOptions;
    if(numExtraOptions > 1) {
        printf("-E- Too many parameters on the command line\n");
        exit(1);
    }

    // print usage if no args given
    if(argc==1) {
        clo_printUsage(list);
        exit(0);
    }

    // list all products
    if(clo_getBool(list, "-l")) {
        printProductList();
        exit(0);
    }

    // get the sensor ID
    const char* sensorName = clo_getString(list, "sensor");
    int sensorId = sensorName2SensorId(sensorName);
    if(sensorId == -1) {
        printf("-E- Could not find sensor \"%s\"\n", sensorName);
        exit(1);
    }
    // just in case a sensor ID was passed in
    sensorName = sensorId2SensorName(sensorId);

    // list all products recursively
    if(clo_getBool(list, "-r")) {
        printProductListSensor(sensorId);
        exit(0);
    }

    // list detailed info for a single product
    if(numExtraOptions == 0) {
        printf("-E- A product name needs to be given on the command line\n");
        exit(1);
    }

    clo_option_t* option = clo_getOption(list, initialNumOptions);
    char* productName = option->key;

    productInfo_t* productInfo = allocateProductInfo();
    if(findProductInfo(productName, sensorId, productInfo)) {
        printf("sensorName=%s\n", sensorName);
        printProductInfo(productName, productInfo);
        return 0;
    }

    printf("-E- Product \"%s\" is not a valid product\n", productName);
    return 1;
}
