# $1 parameter allows to specify a prefix to identify your tmp folders
function get_tmp_dir(){
  local tmp_dir; tmp_dir=$(mktemp -d -t "$1"_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXX)
  echo $tmp_dir
}

# in debug mode ($1=1), do not delete the temporary directory passed as $2
function clean_tmp_dir(){
  if (( $1==0 )); then
    rm -rf "$2"
  fi
}