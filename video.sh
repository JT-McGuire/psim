#! /bin/bash
#ffmpeg -video_size 5388x2100 -framerate 30 -f x11grab -i :0.0+150,460 -c:v hevc_nvenc -preset lossless -pixel_format yuv420p recording.mp4

# ffmpeg -video_size 3840x2160 -framerate 60 \
#     -hwaccel cuda -hwaccel_output_format cuda \
#     -f x11grab  -i :1.0+0,0 -c:v hevc_nvenc -preset lossless -rc:v vbr_hq -b:v 0 \
#     -pixel_format cuda -threads 12 -vf 'format=nv12,hwupload_cuda' ~/Video/screencap/broll0.mkv

ffmpeg -video_size 3680x2100 -framerate 30 \
    -f x11grab -i :1.0+160,560 \
    -vf 'format=nv12,hwupload_cuda' \
    -c:v hevc_nvenc \
    -preset lossless \
    ~/Video/screencap/broll0.mkv

# ffmpeg -video_size 3680x2100 -framerate 60 \
#     -f x11grab -i :1.0+160,560 \
#     -vf 'format=nv12,hwupload_cuda' \
#     -c:v hevc_nvenc \
#     -preset p1 \
#     -rc:v vbr \
#     -cq:v 20 \
#     -b:v 50M \
#     -bf 0 \
#     -g 60 \
#     -rc-lookahead 8 \
#     -threads 0 \
#     ~/Video/screencap/broll0.mkv





