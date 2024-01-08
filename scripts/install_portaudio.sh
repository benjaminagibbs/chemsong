#!/bin/bash
{
   sudo yum install -y gcc-c++ alsa-lib-devel
   wget http://www.portaudio.com/archives/pa_stable_v190700_20210406.tgz
   tar -xzf pa_stable_v190700_20210406.tgz
   cd portaudio
   ./configure && make
   sudo make install
} &> /var/log/install_portaudio.log