/* Caleb R Tulloss
   Brown Neuromotion Lab
   Eye Orientation Tracker
   Analog Coupling Detection V2 and Motor Control
   Software Version 3
   3/18/2018
*/

#include <SPI.h>
#include <Stepper.h>
#include <math.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BME280.h>

/****************************************************************************/
// STEPPERS

# define STEPS 200

Stepper stepperPhi = Stepper(STEPS, 19,18,17,16);
Stepper stepperTheta = Stepper(STEPS, 15,14,2,3);
/****************************************************************************/
// WEATHER STATION

Adafruit_BME280 bme; // I2C

float temp;
float pressure;
float humidity;

/****************************************************************************/
// PINS

// analog inputs from detection circuit
const uint8_t ch0 = 0;
const uint8_t ch1 = 1;
const uint8_t ch2 = 2;
const uint8_t ch3 = 3;

// accelerometer outputs
const uint8_t selectAccel = 52;
const uint8_t accelInterrupt = 53;

/****************************************************************************/
// OTHER CONSTANTS

// Motor control
const uint8_t thetaStepSize = 5;             // 1 step is 0.018 degrees for rods,
const uint8_t phiStepSize = 5;              // roughly 0.0055 degrees for marble
const uint8_t thetaSpeed = 30;                // rpm
const uint8_t phiSpeed = 60;

// Measurement
const uint8_t numSamples = 32;

// Orientation measurement
const uint16_t numAccelSamples = 64;
const uint8_t accelSamplePower = 6;

// Sweep parameters
const uint8_t checkInterval = 20;
const uint8_t thetaNumSteps = 50;             // single-sided
const float phiLowerBound = -0.262;           // -pi/12 = -15 degrees
const float phiUpperBound = 0.628;            // pi/5 = 36 degrees
const float phiCenter = 0.05;

/****************************************************************************/
// GLOBAL VARIABLES

// Data values from final input
uint16_t data0;
uint16_t data1;
uint16_t data2;
uint16_t data3;

uint32_t sum0;
uint32_t sum1;
uint32_t sum2;
uint32_t sum3;

// Motor control
int8_t thetaDirection;
int8_t phiDirection;
uint8_t thetaCount = 0;
int16_t phiCount = 0;

// Orientation
uint8_t accelXLow, accelYLow, accelZLow;
int8_t accelXHigh, accelYHigh, accelZHigh;
int16_t accelDataX, accelDataY, accelDataZ;
int32_t accelSumX, accelSumY, accelSumZ;
int16_t offsetX, offsetY, offsetZ;
float unitX, unitY, unitZ;
float currentPhi;
bool readyForThetaStep;

// Measurement status
bool notDone;

// Time
unsigned long t;

/****************************************************************************/
// HELPER METHODS

// logData: prints the orientation and coupling data, timestamp, and weather sensor data
void logData() 
{
  Serial.print(thetaCount);
  Serial.print(";");
  Serial.print(phiDirection * phiCount);
  Serial.print(";");
  Serial.print(accelSumX);
  Serial.print(";");
  Serial.print(accelSumY);
  Serial.print(";");
  Serial.print(accelSumZ);
  Serial.print(";");
  Serial.print(sum0);
  Serial.print(";");
  Serial.print(sum1);
  Serial.print(";");
  Serial.print(sum2);
  Serial.print(";");
  Serial.print(sum3);
  Serial.print(";");
  Serial.print(t);
  Serial.print(";");
  Serial.print(temp, 4);
  Serial.print(";");
  Serial.print(pressure, 4);
  Serial.print(";");
  Serial.print(humidity, 4);
  Serial.println(";");
}

// measureWeather: measures the current temperature, pressure, and humidity from the Bosch sensor / Adafruit breakout
void measureWeather()
{
  temp = bme.readTemperature();
  pressure = bme.readPressure();
  humidity = bme.readHumidity();
}

// measureAccelOrientation: measures the current X/Y/Z orientation
void measureAccelOrientation()
{
  // First, setup measurement mode, ultralow noise
  digitalWrite(selectAccel, LOW);
  SPI.transfer(0x0A);             // write instruction
  SPI.transfer(0x2D);             // POWER_CTL register
  SPI.transfer(0x22);             // ultralow noise, measurement mode
  digitalWrite(selectAccel, HIGH);
  delay(10);

  accelSumX = 0;
  accelSumY = 0;
  accelSumZ = 0;

  // Then, read the data
  for (int i = 0; i < numAccelSamples; i++)
  {
    while (digitalRead(accelInterrupt) == 0);   // wait for the data to be ready
    
    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x0E);
    accelXLow = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);
    
    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x0F);
    accelXHigh = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);

    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x10);
    accelYLow = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);
    
    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x11);
    accelYHigh = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);

    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x12);
    accelZLow = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);
    
    digitalWrite(selectAccel, LOW);
    SPI.transfer(0x0B);
    SPI.transfer(0x13);
    accelZHigh = SPI.transfer(0x00);
    digitalWrite(selectAccel, HIGH);

    // Combine the two data bytes
    accelDataX = ((int16_t)(accelXHigh) << 8) + (int16_t)accelXLow;
    accelDataY = ((int16_t)(accelYHigh) << 8) + (int16_t)accelYLow;
    accelDataZ = ((int16_t)(accelZHigh) << 8) + (int16_t)accelZLow;

    // Subtract offset
    accelDataX -= offsetX;
    accelDataY -= offsetY;
    accelDataZ -= offsetZ;

    // Add to total for averaging
    accelSumX += accelDataX;
    accelSumY += accelDataY;
    accelSumZ += accelDataZ;
  }

  digitalWrite(selectAccel, LOW);
  SPI.transfer(0x0A);             // write instruction
  SPI.transfer(0x2D);             // POWER_CTL register
  SPI.transfer(0x00);             // standby
  digitalWrite(selectAccel, HIGH);
  delay(10);
}

// for converting accelerometer data to actual orientation
// standard unit vector calculation
void unitOrientation (int32_t x, int32_t y, int32_t z)
{
  float xFl = (float)x;
  float yFl = (float)y;
  float zFl = (float)z;

  float total = sqrt(sq(xFl) + sq(yFl) + sq(zFl));

  unitX = xFl / total;
  unitY = yFl / total;
  unitZ = zFl / total;
}

/****************************************************************************/
// MAIN METHODS

void setup()
{
  // Pin modes
  pinMode(selectAccel, OUTPUT);
  pinMode(accelInterrupt, INPUT);

  // Initialize SPI
  SPI.begin();
  SPI.setDataMode(SPI_MODE0); //CPHA = CPOL = 0    MODE = 0
  delay(1000);

  // Accelerometer soft reset
  digitalWrite(selectAccel, LOW);
  SPI.transfer(0x0A);             // write instruction
  SPI.transfer(0x1F);             // SOFT_RESET register
  SPI.transfer(0x52);             // letter "R"
  digitalWrite(selectAccel, HIGH);
  delay(2000);                    // need to wait at least 500

  digitalWrite(selectAccel, LOW);
  SPI.transfer(0x0A);             // write instruction
  SPI.transfer(0x2A);             // INTMAP1 register
  SPI.transfer(0x01);             // setup data ready interrupt
  digitalWrite(selectAccel, HIGH);
  
  // Set ADC and PWM to be 12-bit
  analogReadResolution(12);
  analogWriteResolution(12);
  
  // Begin serial communication
  Serial.begin(115200);

  // Motor setup
  stepperTheta.setSpeed(thetaSpeed);
  stepperPhi.setSpeed(phiSpeed);
  thetaDirection = -1;
  phiDirection = 1;

  Serial.println("SETUP COMPLETE");

  // Calibrate accelerometer
  offsetX = 0;
  offsetY = 0;
  offsetZ = 0;
  measureAccelOrientation();
  Serial.print("SumX: ");
  Serial.print(accelSumX);
  Serial.print("\tSumY: ");
  Serial.print(accelSumY);
  Serial.print("\tSumZ: ");
  Serial.println(accelSumZ);
  // average
  offsetX = (int16_t)(accelSumX >> accelSamplePower);
  offsetY = (int16_t)(accelSumY >> accelSamplePower);
  offsetZ = (int16_t)(accelSumZ >> accelSamplePower) + 1024;
  // +1024 because in calibration position, Z acceleration should be -1g
  Serial.print("OffsetX: ");
  Serial.print(offsetX);
  Serial.print("\tOffsetY: ");
  Serial.print(offsetY);
  Serial.print("\tOffsetZ: ");
  Serial.println(offsetZ);

  // Setup weather station
  bool status;
  status = bme.begin();
  if (!status) {
    Serial.println("Could not find a valid BME280 sensor, check wiring!");
    while (1);
  }

  // We are not done
  notDone = true;
}

void loop()
{
  if (notDone)
  {
    delay(200);
    sum0 = 0;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    // Measure amplified data values
    for (uint8_t i = 0; i < numSamples; i++)
    {
      
      data0 = analogRead(ch0);
      data1 = analogRead(ch1);
      data2 = analogRead(ch2);
      data3 = analogRead(ch3);

      sum0 += data0;
      sum1 += data1;
      sum2 += data2;
      sum3 += data3;
    }

    // measure ground truth and weather data
    measureAccelOrientation();
    measureWeather();
    // timestamp for loggging the data
    t = millis();
    // Store the data
    logData();

    // Turn eye to new orientation
    if (thetaCount == thetaNumSteps)  // if at end of theta range, we're done
    {
      Serial.println("done!");
      notDone = false;
    }
    else
    {
      if (phiCount == checkInterval)  // every <checkInterval> steps, check where we are
      {
        // phi is the determining factor
        unitOrientation(accelSumX, accelSumY, accelSumZ);
        currentPhi = asin(unitX);
        phiCount = 0;
        
        // if we are in the center and moving in negative direction
        if (abs(currentPhi) < phiCenter && readyForThetaStep)         
        {
          delay(1000);
          stepperTheta.step(thetaDirection * thetaStepSize);    // increment theta
          thetaCount++;
          delay(1000);
          readyForThetaStep = false;
        }
        
        else
        {
          // if phi is past its upper bound
          if (currentPhi > phiUpperBound)
          {
            phiDirection = 1;           // switch direction
            readyForThetaStep = true;
          }
          // if phi is past its lower bound
          else if (currentPhi < phiLowerBound)
          {
            phiDirection = -1;          // switch direction
            readyForThetaStep = false;
          }
        }
      }
      
      // step in the phi direction
      stepperPhi.step(phiDirection * phiStepSize);
      phiCount++;
    }
  }
}

