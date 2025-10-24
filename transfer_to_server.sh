rsync -avz --progress --checksum \
  --exclude 'build/' \
  --exclude '.ccache/' \
  --exclude '*.log' \
  --exclude 'z_output_data_/' \
  --exclude 'third-party' \
  --itemize-changes \
  mutable/ oczerwia@so014:/home/oczerwia/mutable/
