while [ $# -gt 0 ]; do

pdf=$1
echo "Converting "$pdf" ..."
pngfile=`echo "$pdf" | sed 's/\.\w*$/.png/'`
inkscape "$pdf" -z --export-dpi=600 --export-area-drawing --export-png="$pngfile"
echo "Converted to "$pngfile""
shift

done

echo "All jobs done. Exiting."
