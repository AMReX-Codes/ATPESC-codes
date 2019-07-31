#/usr/bin/env bash

# use ffmpeg to convert an mp4 to a gif
# run as: ./mp4_to_gif.sh [file name w/o extension] [framerate]

ffmpeg -i $1.mp4 -r $2 $1.gif
