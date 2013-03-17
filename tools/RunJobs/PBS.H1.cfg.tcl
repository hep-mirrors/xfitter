array set Farm {
  Type "PBS"
  Name "H1"
  Queue "farm12H1G"
  Driver jpbs.sh
  TimeLimit "0:15:00"
  WorkDir /pbstmp
}
