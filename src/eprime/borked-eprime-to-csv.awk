#!/usr/bin/awk

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

/Age:/ { age=$2 ; block=1 }
/Sex:/ { sex=$2 }
/Handedness:/ { hand=$2 }
/Gender:/ { face_gender=$2 }
/Emotion:/ { emotion=$2 }
/GazeDirection:/ { gaze=$2 }
/TargetPosition:/ { target=$2 }
/Target.RTTime:/ { button_time=$2 }

/Target.OnsetTime:/ {
    # extrapolate missing FixationPoint.OnsetTime
    print $2-2000 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n"
    # extrapolate missing DirectGaze.OnsetTime
    print $2-1200 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" "n" "\t" "n" "\t" "n"
    # extrapolate missing AvertedGaze.OnsetTime
    print $2-200 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" "n" "\t" "n"

    print $2 "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" target "\t" "n"
    onset_time=$2
}

/Target.RESP:/ {
    if (button_time != 0) { # don't assume a button was actually pressed
        resp = (($2 == "a" || $2 == "b") ? "left" : "right")
        print button_time "\t" age "\t" sex "\t" hand "\t" block "\t" "obj" "\t" "face" "\t" face_gender "\t" emotion "\t" gaze "\t" target "\t" resp
    }

    # extrapolate missing FinishTiral.OnsetTime
    print onset_time+1500 "\t" age "\t" sex "\t" hand "\t" block "\t" "blank" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n" "\t" "n"
    block=block+1
}

