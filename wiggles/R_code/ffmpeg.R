ffmpeg <- function(folder="~/spiral_pngs", file="spirals.mp4", framerate=30, firstframe=1, lastframe=5) {
  input_pattern <- path.expand(file.path(folder, "*.png"))
  output_file <- path.expand(file.path(folder, file))
  
  command <- paste0("ffmpeg -framerate ", framerate, 
    " -pattern_type glob -y -i \"", input_pattern, # the folder with pngs to be made into a video
    "\" -vf \"tpad=start_mode=clone:start_duration=", firstframe, # enable a pause at the start... 
    ":stop_mode=clone:stop_duration=", lastframe, # and at the end
    ",scale=trunc(iw/2)*2:trunc(ih/2)*2\" ", # avoid problems with frame dimensions
    "-c:v libx264 -profile:v high -level:v 5.1 -crf 18 -preset slow -pix_fmt yuv420p ", # most software should be able to play these videos. -level 4.0 would be safer though for older devices. -crf 18 gives high-quality at the cost of large file size. -preset slow to help retain more detail upon compressin
    "-movflags +faststart -an \"", output_file, "\"") # optimise for web playing. No audio
system(command)
}



combine.av <- function(folder, video, audio, output) {
  video_fl <- path.expand(file.path(folder, video))
  audio_fl <- path.expand(file.path(folder, audio))	
  output_fl <- path.expand(file.path(folder, output))
  
  command <- paste0("ffmpeg -y -i ", video_fl, " -i ", 
    audio_fl, " -c:v copy -ac 1 -c:a aac -shortest ", output_fl)
  system(command)		
}

