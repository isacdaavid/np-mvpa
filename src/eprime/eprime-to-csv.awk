#!/usr/bin/awk

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

# build pyMVPA-friendly design matrix from e-prime event list.

# CSV output format:
#
# the matrix is constructed by searching for all the occurences of certain
# e-prime variables (here preppended with '$'). variables listed on the first
# column trigger the creation of a new row, which is filled according to the
# current state of the program and some default values.
#
# TIME                     AGE  SEX  HANDEDNESS  BLOCKNUM VISUA FACE FACE_GENDER EMOTION  GAZE           TARGET          RESPONSE
#
# $FixationPoint.OnsetTime $Age $Sex $Handedness #        obj   n    n           n        n              n               n
# $DirectGaze.OnsetTime    $Age $Sex $Handedness #        obj   y    $Gender     $Emotion n              n               n
# $AvertedGaze.OnsetTime   $Age $Sex $Handedness #        obj   y    $Gender     $Emotion $GazeDirection n               n
# $Target.OnsetTime        $Age $Sex $Handedness #        obj   y    $Gender     $Emotion $GazeDirection $TargetPosition n
# $Target.RTTime           $Age $Sex $Handedness #        obj   y    $Gender     $Emotion $GazeDirection $TargetPosition $Target.RESP
# $FinishTiral.OnsetTime   $Age $Sex $Handedness #        blank n    n           n        n              n               n

/Age:/ { age=$2 ; block=1 ; bad_file=1 }
/Sex:/ { sex=$2 }
/Handedness:/ { hand=$2 }
/Gender:/ { face_gender=$2 }
/Emotion:/ { emotion=$2 }
/GazeDirection:/ { gaze=$2 }
/TargetPosition:/ { target=$2 }
/Target.RTTime:/ { button_time=$2 }

/FixationPoint.OnsetTime:/ {
    bad_file=0
    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n"
}

/DirectGaze.OnsetTime:/ {
    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" "n" "\t" "n" "\t" "n"
}

/AvertedGaze.OnsetTime:/ {
    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" "n" "\t" "n"
}

/Target.OnsetTime:/ {
    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" target "\t" "n"
}

/Target.RESP:/ {
    if (button_time != 0) { # don't assume a button was actually pressed
        resp = (($2 == "a" || $2 == "b") ? "left" : "right")
        print button_time "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" target "\t" resp
    }
    if (bad_file == 1) {
	block=block+1
    }
}

/FinishTiral.OnsetTime:/ {
    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "blank" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n"
    block=block+1
}

# direct gaze sequence (gaze cuing control) ####################################

# onset times are borked, commented out for now

# /FixationTr.OnsetTime:/ {
#     print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n"
# }

# /DirectGazeNull.OnsetTime:/ {
#     print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" "n" "\t" "n" "\t" "n"
# }
