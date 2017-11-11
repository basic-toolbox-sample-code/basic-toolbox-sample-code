
echo Start profiling
 perf record -a -e cycles:pp -- sleep 10000 &

#perf stat -a  -- sleep 10000 &

#amplxe-cl -collect advanced-hotspots -knob sampling-interval=1 -knob collection-detail=hotspots-sampling -knob event-mode=all -knob enable-user-tasks=false -knob enable-gpu-usage=false -knob enable-gpu-runtimes=false -target-duration-type=veryshort -data-limit=500 -slow-frames-threshold=40 -fast-frames-threshold=100 -no-analyze-kvm-guest --duration unlimited &


#amplxe-cl -collect hotspots -knob sampling-interval=10 -knob enable-user-tasks=true -knob enable-gpu-usage=false -knob enable-gpu-runtimes=false -follow-child -mrte-mode=auto -target-duration-type=veryshort -no-analyze-system -data-limit=500 -slow-frames-threshold=40 -fast-frames-threshold=100 -no-analyze-kvm-guest --target-process pbfs.x &

#sleep 5
