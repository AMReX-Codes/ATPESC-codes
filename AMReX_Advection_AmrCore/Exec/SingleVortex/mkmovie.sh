#!/usr/bin/bash
#for f in $(ls -d plt?????); do ./fsnapshot.gnu.ex -v phi $f; convert $f.phi.ppm $f.phi.png; rm $f.phi.ppm; done

for f in $(ls -d plt?????); do
    python3 slice.py $f -f phi -ax z -min 0.0 -max 1.0 &
done

wait

python3 ffmpeg_make_mp4 plt*.png -s -ifps 5 > /dev/null 2>&1
