#!/usr/bin/bash
for f in $(ls -d plt?????); do ./fsnapshot.gnu.ex -v phi $f; convert $f.phi.ppm $f.phi.png; rm $f.phi.ppm; done
python3 ffmpeg_make_mp4 plt*.phi.png -s -ifps 1
