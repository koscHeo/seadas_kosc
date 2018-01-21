/*
*******************************************************************************
* hico_L0_hdr_echo.c
*
* Read the HICO L0 header and echo the values as text
*
* See usage_call at the end of this document or hico.h for the breakdown
* of the header values.
*
* exit codes:
*   0 = normal
*   1 = improper usage
*   2 = file not accessible
*
* Originally written by Karen Patterson
* November 2009
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hico.h"

void usage_call(void);
char *which_endian();
void *swap_endian(void* Addr, const int Nbytes);
void *bit_string_lsb8(unsigned short int intnum);
void *bit_string_lsb16(unsigned int intnum);
unsigned long int bit_int(int exponent);

int main (int argc, char *argv[])
{
  FILE *in_file;
  char *in_filename;
  struct HicoL0Header hico_header;
  int i;

  /* only run if there is the proper number of arguments */
  if (argc == 2)
    {
    in_filename = argv[1];
    }
  else
    usage_call();

  /* only run if the input file is readable */
  if ((in_file = fopen(in_filename, "rb")) != NULL)
    {
    /* read the header */
    fread(&hico_header, sizeof hico_header, 1, in_file);
    /* The header is big-endian.  The data is little-endian. */
    /* if we're on a little-endian machine, byte swap when necessary */
    if (which_endian() == "LITTLE_ENDIAN")
      {
      hico_header.PkgLen = *(unsigned int *) swap_endian((void*)&hico_header.PkgLen, sizeof hico_header.PkgLen);
      hico_header.FFSubSec = *(unsigned int *) swap_endian((void*)&hico_header.FFSubSec, sizeof hico_header.FFSubSec);
      /* need to strip off the first 12 bits - they should be zeroes, but don't rely on it */
      /* so the value should be less than 1048576 */
      for(i=31; i>=20; i--)
        {
        if (hico_header.FFSubSec >= bit_int(i))
          hico_header.FFSubSec = hico_header.FFSubSec - bit_int(i);
        }
      hico_header.FFPPSSec = *(unsigned short int *) swap_endian((void*)&hico_header.FFPPSSec, sizeof hico_header.FFPPSSec);
      hico_header.FFPPSSubSec = *(unsigned short int *) swap_endian((void*)&hico_header.FFPPSSubSec, sizeof hico_header.FFPPSSubSec);
      hico_header.LFPPSSec = *(unsigned short int *) swap_endian((void*)&hico_header.LFPPSSec, sizeof hico_header.LFPPSSec);
      hico_header.LFPPSSubSec = *(unsigned short int *) swap_endian((void*)&hico_header.LFPPSSubSec, sizeof hico_header.LFPPSSubSec);
      hico_header.ROIW = *(unsigned short int *) swap_endian((void*)&hico_header.ROIW, sizeof hico_header.ROIW);
      hico_header.ROIH = *(unsigned short int *) swap_endian((void*)&hico_header.ROIH, sizeof hico_header.ROIH);
      hico_header.ROIX = *(unsigned short int *) swap_endian((void*)&hico_header.ROIX, sizeof hico_header.ROIX);
      hico_header.ROIY = *(unsigned short int *) swap_endian((void*)&hico_header.ROIY, sizeof hico_header.ROIY);
      hico_header.TotalFrames = *(unsigned short int *) swap_endian((void*)&hico_header.TotalFrames, sizeof hico_header.TotalFrames);
      hico_header.PreDarkFrameCount = *(unsigned short int *) swap_endian((void*)&hico_header.PreDarkFrameCount, sizeof hico_header.PreDarkFrameCount);
      hico_header.SceneFrameCount = *(unsigned short int *) swap_endian((void*)&hico_header.SceneFrameCount, sizeof hico_header.SceneFrameCount);
      hico_header.PostDarkFrameCount = *(unsigned short int *) swap_endian((void*)&hico_header.PostDarkFrameCount, sizeof hico_header.PostDarkFrameCount);
      hico_header.CommandID = *(unsigned short int *) swap_endian((void*)&hico_header.CommandID, sizeof hico_header.CommandID);
      hico_header.ExposureTime = *(unsigned int *) swap_endian((void*)&hico_header.ExposureTime, sizeof hico_header.ExposureTime);
      hico_header.FFPositionX = *(float *) swap_endian((void*)&hico_header.FFPositionX, sizeof hico_header.FFPositionX);
      hico_header.FFPositionY = *(float *) swap_endian((void*)&hico_header.FFPositionY, sizeof hico_header.FFPositionY);
      hico_header.FFPositionZ = *(float *) swap_endian((void*)&hico_header.FFPositionZ, sizeof hico_header.FFPositionZ);
      hico_header.FFVelocityX = *(float *) swap_endian((void*)&hico_header.FFVelocityX, sizeof hico_header.FFVelocityX);
      hico_header.FFVelocityY = *(float *) swap_endian((void*)&hico_header.FFVelocityY, sizeof hico_header.FFVelocityY);
      hico_header.FFVelocityZ = *(float *) swap_endian((void*)&hico_header.FFVelocityZ, sizeof hico_header.FFVelocityZ);
      hico_header.FFQuaternionScalar = *(float *) swap_endian((void*)&hico_header.FFQuaternionScalar, sizeof hico_header.FFQuaternionScalar);
      hico_header.FFQuaternionX = *(float *) swap_endian((void*)&hico_header.FFQuaternionX, sizeof hico_header.FFQuaternionX);
      hico_header.FFQuaternionY = *(float *) swap_endian((void*)&hico_header.FFQuaternionY, sizeof hico_header.FFQuaternionY);
      hico_header.FFQuaternionZ = *(float *) swap_endian((void*)&hico_header.FFQuaternionZ, sizeof hico_header.FFQuaternionZ);
      hico_header.FFGNCSec = *(unsigned int *) swap_endian((void*)&hico_header.FFGNCSec, sizeof hico_header.FFGNCSec);
      hico_header.TotalErrCount = *(unsigned short int *) swap_endian((void*)&hico_header.TotalErrCount, sizeof hico_header.TotalErrCount);
      hico_header.LFPositionX = *(float *) swap_endian((void*)&hico_header.LFPositionX, sizeof hico_header.LFPositionX);
      hico_header.LFPositionY = *(float *) swap_endian((void*)&hico_header.LFPositionY, sizeof hico_header.LFPositionY);
      hico_header.LFPositionZ = *(float *) swap_endian((void*)&hico_header.LFPositionZ, sizeof hico_header.LFPositionZ);
      hico_header.LFVelocityX = *(float *) swap_endian((void*)&hico_header.LFVelocityX, sizeof hico_header.LFVelocityX);
      hico_header.LFVelocityY = *(float *) swap_endian((void*)&hico_header.LFVelocityY, sizeof hico_header.LFVelocityY);
      hico_header.LFVelocityZ = *(float *) swap_endian((void*)&hico_header.LFVelocityZ, sizeof hico_header.LFVelocityZ);
      hico_header.LFQuaternionScalar = *(float *) swap_endian((void*)&hico_header.LFQuaternionScalar, sizeof hico_header.LFQuaternionScalar);
      hico_header.LFQuaternionX = *(float *) swap_endian((void*)&hico_header.LFQuaternionX, sizeof hico_header.LFQuaternionX);
      hico_header.LFQuaternionY = *(float *) swap_endian((void*)&hico_header.LFQuaternionY, sizeof hico_header.LFQuaternionY);
      hico_header.LFQuaternionZ = *(float *) swap_endian((void*)&hico_header.LFQuaternionZ, sizeof hico_header.LFQuaternionZ);
      hico_header.LFGNCSec = *(unsigned int *) swap_endian((void*)&hico_header.LFGNCSec, sizeof hico_header.LFGNCSec);
      hico_header.EventFlags = *(unsigned short int *) swap_endian((void*)&hico_header.EventFlags, sizeof hico_header.EventFlags);
      hico_header.PreDarkOffset = *(unsigned int *) swap_endian((void*)&hico_header.PreDarkOffset, sizeof hico_header.PreDarkOffset);
      hico_header.SceneOffset = *(unsigned int *) swap_endian((void*)&hico_header.SceneOffset, sizeof hico_header.SceneOffset);
      hico_header.PostDarkOffset = *(unsigned int *) swap_endian((void*)&hico_header.PostDarkOffset, sizeof hico_header.PostDarkOffset);
      hico_header.AnalogAmpGain = *(unsigned int *) swap_endian((void*)&hico_header.AnalogAmpGain, sizeof hico_header.AnalogAmpGain);
      hico_header.TimeOnInterval = *(unsigned short int *) swap_endian((void*)&hico_header.TimeOnInterval, sizeof hico_header.TimeOnInterval);
      hico_header.TimeDarkInterval = *(unsigned short int *) swap_endian((void*)&hico_header.TimeDarkInterval, sizeof hico_header.TimeDarkInterval);
      hico_header.TimeObsInterval = *(unsigned short int *) swap_endian((void*)&hico_header.TimeObsInterval, sizeof hico_header.TimeObsInterval);
      hico_header.TimeAngleInterval = *(unsigned short int *) swap_endian((void*)&hico_header.TimeAngleInterval, sizeof hico_header.TimeAngleInterval);
      hico_header.TriggerPeriodCount = *(unsigned short int *) swap_endian((void*)&hico_header.TriggerPeriodCount, sizeof hico_header.TriggerPeriodCount);
      hico_header.EMGain = *(unsigned short int *) swap_endian((void*)&hico_header.EMGain, sizeof hico_header.EMGain);
      hico_header.PreDarkTelemetry = *(unsigned short int *) swap_endian((void*)&hico_header.PreDarkTelemetry, sizeof hico_header.PreDarkTelemetry);
      hico_header.PreDarkEncoderPos = *(short int *) swap_endian((void*)&hico_header.PreDarkEncoderPos, sizeof hico_header.PreDarkEncoderPos);
      hico_header.PreDarkTriggerCount = *(unsigned short int *) swap_endian((void*)&hico_header.PreDarkTriggerCount, sizeof hico_header.PreDarkTriggerCount);
      hico_header.SceneTelemetry = *(unsigned short int *) swap_endian((void*)&hico_header.SceneTelemetry, sizeof hico_header.SceneTelemetry);
      hico_header.SceneEncoderPos = *(short int *) swap_endian((void*)&hico_header.SceneEncoderPos, sizeof hico_header.SceneEncoderPos);
      hico_header.SceneTriggerCount = *(unsigned short int *) swap_endian((void*)&hico_header.SceneTriggerCount, sizeof hico_header.SceneTriggerCount);
      hico_header.PostDarkTelemetry = *(unsigned short int *) swap_endian((void*)&hico_header.PostDarkTelemetry, sizeof hico_header.PostDarkTelemetry);
      hico_header.PostDarkEncoderPos = *(short int *) swap_endian((void*)&hico_header.PostDarkEncoderPos, sizeof hico_header.PostDarkEncoderPos);
      hico_header.PostDarkTriggerCount = *(unsigned short int *) swap_endian((void*)&hico_header.PostDarkTriggerCount, sizeof hico_header.PostDarkTriggerCount);
      hico_header.TimeMotor = *(unsigned short int *) swap_endian((void*)&hico_header.TimeMotor, sizeof hico_header.TimeMotor);
      hico_header.Spare100 = *(unsigned short int *) swap_endian((void*)&hico_header.Spare100, sizeof hico_header.Spare100);
      hico_header.Spare101 = *(unsigned short int *) swap_endian((void*)&hico_header.Spare101, sizeof hico_header.Spare101);
      hico_header.Spare102 = *(unsigned short int *) swap_endian((void*)&hico_header.Spare102, sizeof hico_header.Spare102);
      hico_header.ErrLogCount = *(unsigned short int *) swap_endian((void*)&hico_header.ErrLogCount, sizeof hico_header.ErrLogCount);
      hico_header.ErrSW0 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW0, sizeof hico_header.ErrSW0);
      hico_header.ErrID0 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID0, sizeof hico_header.ErrID0);
      hico_header.ErrSW1 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW1, sizeof hico_header.ErrSW1);
      hico_header.ErrID1 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID1, sizeof hico_header.ErrID1);
      hico_header.ErrSW2 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW2, sizeof hico_header.ErrSW2);
      hico_header.ErrID2 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID2, sizeof hico_header.ErrID2);
      hico_header.ErrSW3 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW3, sizeof hico_header.ErrSW3);
      hico_header.ErrID3 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID3, sizeof hico_header.ErrID3);
      hico_header.ErrSW4 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW4, sizeof hico_header.ErrSW4);
      hico_header.ErrID4 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID4, sizeof hico_header.ErrID4);
      hico_header.ErrSW5 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW5, sizeof hico_header.ErrSW5);
      hico_header.ErrID5 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID5, sizeof hico_header.ErrID5);
      hico_header.ErrSW6 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW6, sizeof hico_header.ErrSW6);
      hico_header.ErrID6 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID6, sizeof hico_header.ErrID6);
      hico_header.ErrSW7 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW7, sizeof hico_header.ErrSW7);
      hico_header.ErrID7 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID7, sizeof hico_header.ErrID7);
      hico_header.ErrSW8 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW8, sizeof hico_header.ErrSW8);
      hico_header.ErrID8 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID8, sizeof hico_header.ErrID8);
      hico_header.ErrSW9 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW9, sizeof hico_header.ErrSW9);
      hico_header.ErrID9 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID9, sizeof hico_header.ErrID9);
      hico_header.ErrSW10 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW10, sizeof hico_header.ErrSW10);
      hico_header.ErrID10 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID10, sizeof hico_header.ErrID10);
      hico_header.ErrSW11 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrSW11, sizeof hico_header.ErrSW11);
      hico_header.ErrID11 = *(unsigned short int *) swap_endian((void*)&hico_header.ErrID11, sizeof hico_header.ErrID11);
      }

    /* print the header values to stdout */
    /* in the structure, the Sync word is not a terminated string */
    fprintf(stdout, "Header Sync Word = ");
    for(i=0; i<strlen(hico_header.Sync)-1; i++)
      fprintf(stdout, "%c",hico_header.Sync[i]);
    fprintf(stdout,"\n");
    fprintf(stdout, "Data Package Length = %u\n",hico_header.PkgLen);
    fprintf(stdout, "First Frame Time Preamble = %X\n",hico_header.FFTimePre);
    fprintf(stdout, "First Frame BCD Year MSB = %02X\n",hico_header.FFYearMSB);
    fprintf(stdout, "First Frame BCD Year LSB = %02X\n",hico_header.FFYearLSB);
    fprintf(stdout, "First Frame BCD Month = %02X\n",hico_header.FFMonth);
    fprintf(stdout, "First Frame BCD Day = %02X\n",hico_header.FFDay);
    fprintf(stdout, "First Frame BCD Hours = %02X\n",hico_header.FFHour);
    fprintf(stdout, "First Frame BCD Minutes = %02X\n",hico_header.FFMinute);
    fprintf(stdout, "First Frame BCD Seconds = %02X\n",hico_header.FFSecond);
    fprintf(stdout, "First Frame Time Sub-Seconds = %u\n",hico_header.FFSubSec);
    fprintf(stdout, "First Frame PPS Hardware Seconds Count = %u\n",hico_header.FFPPSSec);
    fprintf(stdout, "First Frame PPS Hardware Sub-Seconds Count = %u\n",hico_header.FFPPSSubSec);
    fprintf(stdout, "Last Frame PPS Hardware Seconds Count = %u\n",hico_header.LFPPSSec);
    fprintf(stdout, "Last Frame PPS Hardware Sub-Seconds Count = %u\n",hico_header.LFPPSSubSec);
    fprintf(stdout, "ROI Width = %u\n",hico_header.ROIW);
    fprintf(stdout, "ROI Height = %u\n",hico_header.ROIH);
    fprintf(stdout, "ROI X = %u\n",hico_header.ROIX);
    fprintf(stdout, "ROI Y = %u\n",hico_header.ROIY);
    fprintf(stdout, "Horizontal Binning = %u\n",hico_header.HorizBin);
    fprintf(stdout, "Vertical Binning = %u\n",hico_header.VertBin);
    fprintf(stdout, "Readout Port = %u\n",hico_header.ReadoutPort);
    fprintf(stdout, "Readout Speed = %u\n",hico_header.ReadoutSpeed);
    fprintf(stdout, "ICU Application Software Version = %u\n",hico_header.ICUAppVer);
    fprintf(stdout, "EM GainClearing Mode = %u\n",hico_header.EMClearingMode);
    fprintf(stdout, "Total Frames = %u\n",hico_header.TotalFrames);
    fprintf(stdout, "Pre-Exposure Dark Frame Count = %u\n",hico_header.PreDarkFrameCount);
    fprintf(stdout, "Scene Data Frame Count = %u\n",hico_header.SceneFrameCount);
    fprintf(stdout, "Post-Exposure Dark Frame Count = %u\n",hico_header.PostDarkFrameCount);
    fprintf(stdout, "Command ID = %u\n",hico_header.CommandID);
    fprintf(stdout, "Exposure Time = %u\n",hico_header.ExposureTime);
    fprintf(stdout, "First Frame ISS Position Vector X Component = %0.7E\n",hico_header.FFPositionX);
    fprintf(stdout, "First Frame ISS Position Vector Y Component = %0.7E\n",hico_header.FFPositionY);
    fprintf(stdout, "First Frame ISS Position Vector Z Component = %0.7E\n",hico_header.FFPositionZ);
    fprintf(stdout, "First Frame ISS Velocity Vector X Component = %0.7E\n",hico_header.FFVelocityX);
    fprintf(stdout, "First Frame ISS Velocity Vector Y Component = %0.7E\n",hico_header.FFVelocityY);
    fprintf(stdout, "First Frame ISS Velocity Vector Z Component = %0.7E\n",hico_header.FFVelocityZ);
    fprintf(stdout, "First Frame ISS Attitude Quaternion Scalar Component = %0.7E\n",hico_header.FFQuaternionScalar);
    fprintf(stdout, "First Frame ISS Attitude Quaternion X Component = %0.7E\n",hico_header.FFQuaternionX);
    fprintf(stdout, "First Frame ISS Attitude Quaternion Y Component = %0.7E\n",hico_header.FFQuaternionY);
    fprintf(stdout, "First Frame ISS Attitude Quaternion Z Component = %0.7E\n",hico_header.FFQuaternionZ);
    fprintf(stdout, "First Frame GNC Time Seconds = %u\n",hico_header.FFGNCSec);
    fprintf(stdout, "First Frame GNC Time Sub-Seconds = %u\n",hico_header.FFGNCSubSec);
    fprintf(stdout, "First Frame GNC Quality Flags = %s\n",bit_string_lsb8(hico_header.FFGNCQualFlags));
    fprintf(stdout, "Total Error Count = %u\n",hico_header.TotalErrCount);
    fprintf(stdout, "Last Frame ISS Position Vector X Component = %0.7E\n",hico_header.LFPositionX);
    fprintf(stdout, "Last Frame ISS Position Vector Y Component = %0.7E\n",hico_header.LFPositionY);
    fprintf(stdout, "Last Frame ISS Position Vector Z Component = %0.7E\n",hico_header.LFPositionZ);
    fprintf(stdout, "Last Frame ISS Velocity Vector X Component = %0.7E\n",hico_header.LFVelocityX);
    fprintf(stdout, "Last Frame ISS Velocity Vector Y Component = %0.7E\n",hico_header.LFVelocityY);
    fprintf(stdout, "Last Frame ISS Velocity Vector Z Component = %0.7E\n",hico_header.LFVelocityZ);
    fprintf(stdout, "Last Frame ISS Attitude Quaternion Scalar Component = %0.7E\n",hico_header.LFQuaternionScalar);
    fprintf(stdout, "Last Frame ISS Attitude Quaternion X Component = %0.7E\n",hico_header.LFQuaternionX);
    fprintf(stdout, "Last Frame ISS Attitude Quaternion Y Component = %0.7E\n",hico_header.LFQuaternionY);
    fprintf(stdout, "Last Frame ISS Attitude Quaternion Z Component = %0.7E\n",hico_header.LFQuaternionZ);
    fprintf(stdout, "Last Frame GNC Time Seconds = %u\n",hico_header.LFGNCSec);
    fprintf(stdout, "Last Frame GNC Time Sub-Seconds = %u\n",hico_header.LFGNCSubSec);
    fprintf(stdout, "Last Frame GNC Quality Flags = %s\n",bit_string_lsb8(hico_header.LFGNCQualFlags));
    fprintf(stdout, "Event Data Set Flags = %s\n",bit_string_lsb16(hico_header.EventFlags));
    fprintf(stdout, "Pre-exposure Dark Frame Offset = %u\n",hico_header.PreDarkOffset);
    fprintf(stdout, "Scene Data Frame Offset = %u\n",hico_header.SceneOffset);
    fprintf(stdout, "Post-exposure Dark Frame Offset = %u\n",hico_header.PostDarkOffset);
    fprintf(stdout, "Analog Amplifier Gain = %u\n",hico_header.AnalogAmpGain);
    fprintf(stdout, "Time On Interval = %u\n",hico_header.TimeOnInterval);
    fprintf(stdout, "Time Dark Interval = %u\n",hico_header.TimeDarkInterval);
    fprintf(stdout, "Time Observation Interval = %u\n",hico_header.TimeObsInterval);
    fprintf(stdout, "Time Angle Interval = %u\n",hico_header.TimeAngleInterval);
    fprintf(stdout, "Trigger Pulse Period Count = %u\n",hico_header.TriggerPeriodCount);
    fprintf(stdout, "EM Gain = %u\n",hico_header.EMGain);
    fprintf(stdout, "Pre-exposure Dark Frame Digital Telemetry = %s\n",bit_string_lsb16(hico_header.PreDarkTelemetry));
    fprintf(stdout, "Pre-exposure Dark Frame Encoder Position = %d\n",hico_header.PreDarkEncoderPos);
    fprintf(stdout, "Pre-exposure Dark Frame Trigger Pulse Count = %u\n",hico_header.PreDarkTriggerCount);
    fprintf(stdout, "Scene Data Frame Digital Telemetry = %s\n",bit_string_lsb16(hico_header.SceneTelemetry));
    fprintf(stdout, "Scene Data Frame Encoder Position = %d\n",hico_header.SceneEncoderPos);
    fprintf(stdout, "Scene Data Frame Trigger Pulse Count = %u\n",hico_header.SceneTriggerCount);
    fprintf(stdout, "Post-exposure Dark Frame Digital Telemetry = %s\n",bit_string_lsb16(hico_header.PostDarkTelemetry));
    fprintf(stdout, "Post-exposure Dark Frame Encoder Position = %d\n",hico_header.PostDarkEncoderPos);
    fprintf(stdout, "Post-exposure Dark Frame Trigger Pulse Count = %u\n",hico_header.PostDarkTriggerCount);
    fprintf(stdout, "Time Motor Initialization Interval = %u\n",hico_header.TimeMotor);
    fprintf(stdout, "Spare100 = %u\n",hico_header.Spare100);
    fprintf(stdout, "Spare101 = %u\n",hico_header.Spare101);
    fprintf(stdout, "Spare102 = %u\n",hico_header.Spare102);
    fprintf(stdout, "Error Log Count = %u\n",hico_header.ErrLogCount);
    fprintf(stdout, "Error Status Word 0 = %u\n",hico_header.ErrSW0);
    fprintf(stdout, "Error Frame ID 0 = %u\n",hico_header.ErrID0);
    fprintf(stdout, "Error Status Word 1 = %u\n",hico_header.ErrSW1);
    fprintf(stdout, "Error Frame ID 1 = %u\n",hico_header.ErrID1);
    fprintf(stdout, "Error Status Word 2 = %u\n",hico_header.ErrSW2);
    fprintf(stdout, "Error Frame ID 2 = %u\n",hico_header.ErrID2);
    fprintf(stdout, "Error Status Word 3 = %u\n",hico_header.ErrSW3);
    fprintf(stdout, "Error Frame ID 3 = %u\n",hico_header.ErrID3);
    fprintf(stdout, "Error Status Word 4 = %u\n",hico_header.ErrSW4);
    fprintf(stdout, "Error Frame ID 4 = %u\n",hico_header.ErrID4);
    fprintf(stdout, "Error Status Word 5 = %u\n",hico_header.ErrSW5);
    fprintf(stdout, "Error Frame ID 5 = %u\n",hico_header.ErrID5);
    fprintf(stdout, "Error Status Word 6 = %u\n",hico_header.ErrSW6);
    fprintf(stdout, "Error Frame ID 6 = %u\n",hico_header.ErrID6);
    fprintf(stdout, "Error Status Word 7 = %u\n",hico_header.ErrSW7);
    fprintf(stdout, "Error Frame ID 7 = %u\n",hico_header.ErrID7);
    fprintf(stdout, "Error Status Word 8 = %u\n",hico_header.ErrSW8);
    fprintf(stdout, "Error Frame ID 8 = %u\n",hico_header.ErrID8);
    fprintf(stdout, "Error Status Word 9 = %u\n",hico_header.ErrSW9);
    fprintf(stdout, "Error Frame ID 9 = %u\n",hico_header.ErrID9);
    fprintf(stdout, "Error Status Word 10 = %u\n",hico_header.ErrSW10);
    fprintf(stdout, "Error Frame ID 10 = %u\n",hico_header.ErrID10);
    fprintf(stdout, "Error Status Word 11 = %u\n",hico_header.ErrSW11);
    fprintf(stdout, "Error Frame ID 11 = %u\n",hico_header.ErrID11);
    fclose(in_file);
    }
  else
    {
    fprintf(stderr,"Failed to open input file: %s\n",in_filename);
    exit(2);
    }
  return(0);
}

/* =============================================================
 *    which_endian
 *    determine byte order for this machine
 * ============================================================= */
char *which_endian()
{
  int    int_lword = 0x04030201;
  char *p;

  /* int should be 4 bytes */
  p = (char *)(&int_lword);

  if (p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4) return "LITTLE_ENDIAN";
  if (p[0]==4 && p[1]==3 && p[2]==2 && p[3]==1) return "BIG_ENDIAN";
  if (p[0]==3 && p[1]==4 && p[2]==1 && p[3]==2) return "PDP_ENDIAN";

  fprintf(stderr,"WARNING: Unknown byte order.\n");

  return "UNDEFINED_BYTE_ORDER";
}

/* =============================================================
 *    swap_endian
 *    swap big/little endian for numeric values
 *    not needed for characters since they are 1 byte each
 * ============================================================= */
void *swap_endian(void* Addr, const int Nbytes)
{
  int i;
  static unsigned char swap_val[32];

  switch (Nbytes)
  {
    case 2:
      break;
    case 4:
      break;
    case 8:
      break;
    case 16:
      break;
    case 32:
      break;
    default:
      /* return as-is */
      return (void*)Addr;
  }

  for(i=Nbytes-1;i>=0;i--)
    swap_val[Nbytes-i-1]=*((unsigned char*)Addr+i);

  return (void*)swap_val;
}

/* =============================================================
 *    bit_string_lsb8
 *    returns a string of bits (0's and 1's for a given int)
 *    LSB order, so the 0 bit is on the far right
 * ============================================================= */
void *bit_string_lsb8(unsigned short int intnum)
{
  static char b[] = "00000000";
  int i, n = intnum;

  for(i=0; i<8; i++, n=n/2)
    if (n%2)b[7-i] = '1';
  return (void *)b;
}

/* =============================================================
 *    bit_string_lsb16
 *    returns a string of bits (0's and 1's for a given int)
 *    LSB order, so the 0 bit is on the far right
 * ============================================================= */
void *bit_string_lsb16(unsigned int intnum)
{
  static char b[] = "0000000000000000";
  int i, n = intnum;

  for(i=0; i<16; i++, n=n/2)
    if (n%2) b[15-i] = '1';
  return (void *)b;
}

/* =============================================================
 *    bit_int
 *    returns 2^exponent
 * ============================================================= */
unsigned long int bit_int(int exponent)
{
  int i;
  unsigned long int retval=1;

  for(i=1; i<=exponent; i++)
    retval = retval*2;

  return retval;
}

/* =============================================================
 *    usage_call()
 * ============================================================= */
void usage_call(void)
{
  printf("\n");
  printf("Usage: hico_L0_hdr_echo <HICO_LO_FILE>\n");
  printf("  HICO_L0_FILE = name of HICO L0 .bil file\n");
  printf("\n");
  printf("The HICO L0 .bil file header is returned as a single space separated\n");
  printf("line with the following entries:\n");
  printf("   1 = Header Sync Word\n");
  printf("   2 = Data Package Length\n");
  printf("   3 = First Frame Time Preamble\n");
  printf("   4 = First Frame BCD Year MSB\n");
  printf("   5 = First Frame BCD Year LSB\n");
  printf("   6 = First Frame BCD Month\n");
  printf("   7 = First Frame BCD Day\n");
  printf("   8 = First Frame BCD Hours\n");
  printf("   9 = First Frame BCD Minutes\n");
  printf("  10 = First Frame BCD Seconds\n");
  printf("  11 = First Frame Time Sub-Seconds\n");
  printf("  12 = First Frame PPS Hardware Seconds Count\n");
  printf("  13 = First Frame PPS Hardware Sub-Seconds Count\n");
  printf("  14 = Last Frame PPS Hardware Seconds Count\n");
  printf("  15 = Last Frame PPS Hardware Sub-Seconds Count\n");
  printf("  16 = ROI Width\n");
  printf("  17 = ROI Height\n");
  printf("  18 = ROI X\n");
  printf("  19 = ROI Y\n");
  printf("  20 = Horizontal Binning\n");
  printf("  21 = Vertical Binning\n");
  printf("  22 = Readout Port\n");
  printf("  23 = Readout Speed\n");
  printf("  24 = ICU Application Software Version\n");
  printf("  25 = EM GainClearing Mode\n");
  printf("  26 = Total Frames\n");
  printf("  27 = Pre-Exposure Dark Frame Count\n");
  printf("  28 = Scene Data Frame Count\n");
  printf("  29 = Post-Exposure Dark Frame Count\n");
  printf("  30 = Command ID\n");
  printf("  31 = Exposure Time\n");
  printf("  32 = First Frame ISS Position Vector X Component\n");
  printf("  33 = First Frame ISS Position Vector Y Component\n");
  printf("  34 = First Frame ISS Position Vector Z Component\n");
  printf("  35 = First Frame ISS Velocity Vector X Component\n");
  printf("  36 = First Frame ISS Velocity Vector Y Component\n");
  printf("  37 = First Frame ISS Velocity Vector Z Component\n");
  printf("  38 = First Frame ISS Attitude Quaternion Scalar Component\n");
  printf("  39 = First Frame ISS Attitude Quaternion X Component\n");
  printf("  40 = First Frame ISS Attitude Quaternion Y Component\n");
  printf("  41 = First Frame ISS Attitude Quaternion Z Component\n");
  printf("  42 = First Frame GNC Time Seconds\n");
  printf("  43 = First Frame GNC Time Sub-Seconds\n");
  printf("  44 = First Frame GNC Quality Flags\n");
  printf("  45 = Total Error Count\n");
  printf("  46 = Last Frame ISS Position Vector X Component\n");
  printf("  47 = Last Frame ISS Position Vector Y Component\n");
  printf("  48 = Last Frame ISS Position Vector Z Component\n");
  printf("  49 = Last Frame ISS Velocity Vector X Component\n");
  printf("  50 = Last Frame ISS Velocity Vector Y Component\n");
  printf("  51 = Last Frame ISS Velocity Vector Z Component\n");
  printf("  52 = Last Frame ISS Attitude Quaternion Scalar Component\n");
  printf("  53 = Last Frame ISS Attitude Quaternion X Component\n");
  printf("  54 = Last Frame ISS Attitude Quaternion Y Component\n");
  printf("  55 = Last Frame ISS Attitude Quaternion Z Component\n");
  printf("  56 = Last Frame GNC Time Seconds\n");
  printf("  57 = Last Frame GNC Time Sub-Seconds\n");
  printf("  58 = Last Frame GNC Quality Flags\n");
  printf("  59 = Event Data Set Flags\n");
  printf("  60 = Pre-exposure Dark Frame Offset\n");
  printf("  61 = Scene Data Frame Offset\n");
  printf("  62 = Post-exposure Dark Frame Offset\n");
  printf("  63 = Analog Amplifier Gain\n");
  printf("  64 = Time On Interval\n");
  printf("  65 = Time Dark Interval\n");
  printf("  66 = Time Observation Interval\n");
  printf("  67 = Time Angle Interval\n");
  printf("  68 = Trigger Pulse Period Count\n");
  printf("  69 = EM Gain\n");
  printf("  70 = Pre-exposure Dark Frame Digital Telemetry\n");
  printf("  71 = Pre-exposure Dark Frame Encoder Position\n");
  printf("  72 = Pre-exposure Dark Frame Trigger Pulse Count\n");
  printf("  73 = Scene Data Frame Digital Telemetry\n");
  printf("  74 = Scene Data Frame Encoder Position\n");
  printf("  75 = Scene Data Frame Trigger Pulse Count\n");
  printf("  76 = Post-exposure Dark Frame Digital Telemetry\n");
  printf("  77 = Post-exposure Dark Frame Encoder Position\n");
  printf("  78 = Post-exposure Dark Frame Trigger Pulse Count\n");
  printf("  79 = Time Motor Initialization Interval\n");
  printf("  80 = Spare\n");
  printf("  81 = Spare\n");
  printf("  82 = Spare\n");
  printf("  83 = Error Log Count\n");
  printf("  84 = Error Status Word 0\n");
  printf("  85 = Error Frame ID 0\n");
  printf("  86 = Error Status Word 1\n");
  printf("  87 = Error Frame ID 1\n");
  printf("  88 = Error Status Word 2\n");
  printf("  89 = Error Frame ID 2\n");
  printf("  90 = Error Status Word 3\n");
  printf("  91 = Error Frame ID 3\n");
  printf("  92 = Error Status Word 4\n");
  printf("  93 = Error Frame ID 4\n");
  printf("  94 = Error Status Word 5\n");
  printf("  95 = Error Frame ID 5\n");
  printf("  96 = Error Status Word 6\n");
  printf("  97 = Error Frame ID 6\n");
  printf("  98 = Error Status Word 7\n");
  printf("  99 = Error Frame ID 7\n");
  printf(" 100 = Error Status Word 8\n");
  printf(" 101 = Error Frame ID 8\n");
  printf(" 102 = Error Status Word 9\n");
  printf(" 103 = Error Frame ID 9\n");
  printf(" 104 = Error Status Word 10\n");
  printf(" 105 = Error Frame ID 10\n");
  printf(" 106 = Error Status Word 11\n");
  printf(" 107 = Error Frame ID 11\n");
  exit (1);
}
