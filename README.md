## Setup

Need MATLAB 2011a or above.


## Dependency

You will need a USRP or similar devices to record data in IQ format.

I 	- 16 Bits
Q	- 16 Bits


## Description

Dump 1090 is a Mode S decoder specifically for processing data offline.

The main features are:
- Single bit errors correction using the 24 bit CRC.
- Ability to decode positions and velocities of airplanes
- Decode raw IQ samples from file

## Usage

Run the `dump1090.m`

Specify the parameters as follows:

	seconds = 2;
	dataPath = 'C:\USRP_DATA\ADSB';
	filename = 'adsb_set_10MHz.dat';
	sampRate = 1e7;

The `sampRate` depends on the sampling rate you are using for recording the raw IQ data.
	
You will get output like this:

	 -- We find a packet ... 
	 -- Start index is 18385412 ... 
	 -- DF type is 11 ... 
	 -- Plane ID is 750257 ... 
	 
              ew_dir: 0
         ew_velocity: 136
              ns_dir: 1
         ns_velocity: 258
    vert_rate_source: 0
      vert_rate_sign: 1
           vert_rate: 23
             cprtime: 19292471
            accuracy: 'Error < 0.3m/s'
            velocity: [1x3 double]
               plane: '7502AD'
 

## Credits

Most of the implementation comes from antirez's dump1090: https://github.com/antirez/dump1090/

You can find more powerful implementation at https://github.com/MalcolmRobb/dump1090

## License

Under the MIT license