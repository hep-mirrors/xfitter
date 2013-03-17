array set Farm {
  Type "PBS"
  Name "PL-Grid"
  Queue "plgrid"
  Driver jpbs.sh
  TimeLimit "0:08:00"
  WorkDir /mnt/lustre/scratch/people/<Farm(User)>
}
