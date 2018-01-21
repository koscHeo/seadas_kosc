/*
* hico.h
*
* HICO definitions
*
* Originally written by Karen Patterson
* November 2009
*
*/

#include <stdlib.h>
// #include <endian.h>

/* HICO L0 Header (note: big-endian) */
typedef struct HicoL0Header {
  char Sync[4]; /* Header Sync Word (32 bits) "HICO" */
  unsigned int PkgLen; /* Data Package Length (32 bits) # bytes in file minus 8 */
  unsigned char FFTimePre; /* First Frame Time Preamble (8 bits) hex "50" */
  unsigned char FFYearMSB; /* First Frame BCD Year MSB (8 bits) hex [19-20] */
  unsigned char FFYearLSB; /* First Frame BCD Year LSB (8 bits) hex [00-99] */
  unsigned char FFMonth; /* First Frame BCD Month (8 bits) hex [1-12] */
  unsigned char FFDay; /* First Frame BCD Day (8 bits) hex [1-31] */
  unsigned char FFHour; /* First Frame BCD Hours (8 bits) hex [0-23] */
  unsigned char FFMinute; /* First Frame BCD Minutes (8 bits) hex [0-59] */
  unsigned char FFSecond; /* First Frame BCD Seconds (8 bits) hex [0-59] */
  /* First Frame Time BCD Sub-Seconds (20 bits) hex [0-999,999] */
  /* Read 32 bits and discard the first 12 */
  unsigned int FFSubSec;
  /* some routines read this as (12-bit + 4-bit) 16-bit MSB, 16-bit LSB
  unsigned char Spare8;
  unsigned char FFSubSecMSB;
  unsigned short int FFSubSecLSB;
  */
  unsigned short int FFPPSSec; /* First Frame PPS Hardware Seconds Count (16 bits) */
  unsigned short int FFPPSSubSec; /* First Frame PPS Hardware Sub-Seconds Count (16 bits) */
  unsigned short int LFPPSSec; /* Last Frame PPS Hardware Seconds Count (16 bits) */
  unsigned short int LFPPSSubSec; /* Last Frame PPS Hardware Sub-Seconds Count (16 bits) */
  unsigned short int ROIW; /* ROI Width (16 bits) */
  unsigned short int ROIH; /* ROI Height (16 bits) */
  unsigned short int ROIX; /* ROI X (16 bits) */
  unsigned short int ROIY; /* ROI Y (16 bits) */
  unsigned char HorizBin; /* Horizontal Binning (8 bits) */
  unsigned char VertBin; /* Vertical Binning (8 bits) */
  unsigned char ReadoutPort; /* Readout Port (8 bits) */
  unsigned char ReadoutSpeed; /* Readout Speed (8 bits) */
  unsigned char ICUAppVer; /* ICU Application Software Version (8 bits) */
  unsigned char EMClearingMode; /* EM GainClearing Mode (8 bits) */
  unsigned short int TotalFrames; /* Total Frames (16 bits) */
  unsigned short int PreDarkFrameCount; /* Pre-Exposure Dark Frame Count (16 bits) */
  unsigned short int SceneFrameCount; /* Scene Data Frame Count (16 bits) */
  unsigned short int PostDarkFrameCount; /* Post-Exposure Dark Frame Count (16 bits) */
  unsigned short int CommandID; /* Command ID (16 bits) */
  unsigned int ExposureTime; /* Exposure Time (32 bits) */
  float FFPositionX; /* First Frame ISS Position Vector X Component (32 bits) */
  float FFPositionY; /* First Frame ISS Position Vector Y Component (32 bits) */
  float FFPositionZ; /* First Frame ISS Position Vector Z Component (32 bits) */
  float FFVelocityX; /* First Frame ISS Velocity Vector X Component (32 bits) */
  float FFVelocityY; /* First Frame ISS Velocity Vector Y Component (32 bits) */
  float FFVelocityZ; /* First Frame ISS Velocity Vector Z Component (32 bits) */
  float FFQuaternionScalar; /* First Frame ISS Attitude Quaternion Scalar Component (32 bits) */
  float FFQuaternionX; /* First Frame ISS Attitude Quaternion X Component (32 bits) */
  float FFQuaternionY; /* First Frame ISS Attitude Quaternion Y Component (32 bits) */
  float FFQuaternionZ; /* First Frame ISS Attitude Quaternion Z Component (32 bits) */
  unsigned int FFGNCSec; /* First Frame GNC Time Seconds (32 bits) */
  unsigned char FFGNCSubSec; /* First Frame GNC Time Sub-Seconds (8 bits) */
  unsigned char FFGNCQualFlags; /* First Frame GNC Quality Flags (8 bits) */
  unsigned short int TotalErrCount; /* Total Error Count (16 bits) */
  float LFPositionX; /* Last Frame ISS Position Vector X Component (32 bits) */
  float LFPositionY; /* Last Frame ISS Position Vector Y Component (32 bits) */
  float LFPositionZ; /* Last Frame ISS Position Vector Z Component (32 bits) */
  float LFVelocityX; /* Last Frame ISS Velocity Vector X Component (32 bits) */
  float LFVelocityY; /* Last Frame ISS Velocity Vector Y Component (32 bits) */
  float LFVelocityZ; /* Last Frame ISS Velocity Vector Z Component (32 bits) */
  float LFQuaternionScalar; /* Last Frame ISS Attitude Quaternion Scalar Component (32 bits) */
  float LFQuaternionX; /* Last Frame ISS Attitude Quaternion X Component (32 bits) */
  float LFQuaternionY; /* Last Frame ISS Attitude Quaternion Y Component (32 bits) */
  float LFQuaternionZ; /* Last Frame ISS Attitude Quaternion Z Component (32 bits) */
  unsigned int LFGNCSec; /* Last Frame GNC Time Seconds (32 bits) */
  unsigned char LFGNCSubSec; /* Last Frame GNC Time Sub-Seconds (8 bits) */
  unsigned char LFGNCQualFlags; /* Last Frame GNC Quality Flags (8 bits) */
  unsigned short int EventFlags; /* Event Data Set Flags (16 bits) */
  unsigned int PreDarkOffset; /* Pre-exposure Dark Frame Offset (32 bits) */
  unsigned int SceneOffset; /* Scene Data Frame Offset (32 bits) */
  unsigned int PostDarkOffset; /* Post-exposure Dark Frame Offset (32 bits) */
  unsigned int AnalogAmpGain; /* Analog Amplifier Gain (32 bits) */
  unsigned short int TimeOnInterval; /* Time On Interval (16 bits) */
  unsigned short int TimeDarkInterval; /* Time Dark Interval (16 bits) */
  unsigned short int TimeObsInterval; /* Time Observation Interval (16 bits) */
  unsigned short int TimeAngleInterval; /* Time Angle Interval (16 bits) */
  unsigned short int TriggerPeriodCount; /* Trigger Pulse Period Count (16 bits) */
  unsigned short int EMGain; /* EM Gain (16 bits) */
  unsigned short int PreDarkTelemetry; /* Pre-exposure Dark Frame Digital Telemetry (16 bits) */
  short int PreDarkEncoderPos; /* Pre-exposure Dark Frame Encoder Position (16 bits) */
  unsigned short int PreDarkTriggerCount; /* Pre-exposure Dark Frame Trigger Pulse Count (16 bits) */
  unsigned short int SceneTelemetry; /* Scene Data Frame Digital Telemetry (16 bits) */
  short int SceneEncoderPos; /* Scene Data Frame Encoder Position (16 bits) */
  unsigned short int SceneTriggerCount; /* Scene Data Frame Trigger Pulse Count (16 bits) */
  unsigned short int PostDarkTelemetry; /* Post-exposure Dark Frame Digital Telemetry (16 bits) */
  short int PostDarkEncoderPos; /* Post-exposure Dark Frame Encoder Position (16 bits) */
  unsigned short int PostDarkTriggerCount; /* Post-exposure Dark Frame Trigger Pulse Count (16 bits) */
  unsigned short int TimeMotor; /* Time Motor Initialization Interval (16 bits) */
  unsigned short int Spare100; /* Spare (16 bits) */
  unsigned short int Spare101; /* Spare (16 bits) */
  unsigned short int Spare102; /* Spare (16 bits) */
  unsigned short int ErrLogCount; /* Error Log Count (16 bits) */
  unsigned short int ErrSW0; /* Error Status Word 0 (16 bits) */
  unsigned short int ErrID0; /* Error Frame ID 0 (16 bits) */
  unsigned short int ErrSW1; /* Error Status Word 1 (16 bits) */
  unsigned short int ErrID1; /* Error Frame ID 1 (16 bits) */
  unsigned short int ErrSW2; /* Error Status Word 2 (16 bits) */
  unsigned short int ErrID2; /* Error Frame ID 2 (16 bits) */
  unsigned short int ErrSW3; /* Error Status Word 3 (16 bits) */
  unsigned short int ErrID3; /* Error Frame ID 3 (16 bits) */
  unsigned short int ErrSW4; /* Error Status Word 4 (16 bits) */
  unsigned short int ErrID4; /* Error Frame ID 4 (16 bits) */
  unsigned short int ErrSW5; /* Error Status Word 5 (16 bits) */
  unsigned short int ErrID5; /* Error Frame ID 5 (16 bits) */
  unsigned short int ErrSW6; /* Error Status Word 6 (16 bits) */
  unsigned short int ErrID6; /* Error Frame ID 6 (16 bits) */
  unsigned short int ErrSW7; /* Error Status Word 7 (16 bits) */
  unsigned short int ErrID7; /* Error Frame ID 7 (16 bits) */
  unsigned short int ErrSW8; /* Error Status Word 8 (16 bits) */
  unsigned short int ErrID8; /* Error Frame ID 8 (16 bits) */
  unsigned short int ErrSW9; /* Error Status Word 9 (16 bits) */
  unsigned short int ErrID9; /* Error Frame ID 9 (16 bits) */
  unsigned short int ErrSW10; /* Error Status Word 10 (16 bits) */
  unsigned short int ErrID10; /* Error Frame ID 10 (16 bits) */
  unsigned short int ErrSW11; /* Error Status Word 11 (16 bits) */
  unsigned short int ErrID11; /* Error Frame ID 11 (16 bits) */
} HicoL0Header;
