OpenCL Real-Time Ray Tracer
=================

By Aaron Fritz

If any of my ray tracers are to be truly real time, it's those written with the OpenCL framework. 

This program communicates with an OpenCL-compliant graphics card on a GPGPU-like basis for computing huge amounts of data per second; much more data per second than the central processor can muster. I believe this is attributed to the fact that ray tracing is "embarassingly parallel", and the graphics card may have hundreds or thousands of ~0.5 GHz cores, while the processor has probably four ~3.0 GHz cores. The graphics card is able to process each frame quick enough for it to run at a near-native resolution with a simple scene at sixty frames per second. This was ground-breaking when I first programmed it! OpenCL opened up new doors to having a true real-time ray tracer implementation that I thought would not be viable otherwise.
