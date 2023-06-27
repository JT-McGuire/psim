#! /bin/bash
#ffmpeg -video_size 5388x2100 -framerate 30 -f x11grab -i :0.0+150,460 -c:v hevc_nvenc -preset lossless -pixel_format yuv420p recording.mp4
ffmpeg -video_size 3840x2160 -framerate 60 -hwaccel cuda -hwaccel_output_format cuda -f x11grab -i :0.0+0,910 -c:v hevc_nvenc -preset lossless -pixel_format yuv420p broll0.mp4