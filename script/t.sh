#/bin/bash

del=$1
[ -z $del ] && del=65

cawk -v del=$del '
   BEGIN{ start=systime(); start -= del
          print sec_fmt( systime()-start  )
   }
   function sec_fmt(secs) {
      m = ""; s = ""; h = ""
      if (secs >= 3600) {
         h = int(secs/3600) "h"
         secs -= h * 3600
      }
      if (secs >= 60) {
         m = int(secs/60) "m"
         secs -= m * 60
      }
      if (secs > 0)
         s = secs "s"
      return h m s
   }
'
