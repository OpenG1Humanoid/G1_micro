// G1 micro: the Unitree G1 parts shrunk to 0.35 and edited for micro servos. Keep next to STL/. F5.
// Joint angles in degrees, same names, zero pose and signs as the G1 (elbow 0 = forearm forward).
/* [Waist] */
waist_yaw = 0; // [-150:1:150]
/* [Left] */
left_hip_pitch = 0; // [-145:1:165]
left_hip_roll = 0; // [-30:1:170]
left_hip_yaw = 0; // [-158:1:158]
left_knee = 0; // [-5:1:165]
left_ankle_pitch = 0; // [-50:1:30]
left_ankle_roll = 0; // [-15:1:15]
left_shoulder_pitch = 0; // [-177:1:153]
left_shoulder_roll = 0; // [-91:1:129]
left_shoulder_yaw = 0; // [-150:1:150]
left_elbow = 0; // [-60:1:120]
/* [Right] */
right_hip_pitch = 0; // [-145:1:165]
right_hip_roll = 0; // [-170:1:30]
right_hip_yaw = 0; // [-158:1:158]
right_knee = 0; // [-5:1:165]
right_ankle_pitch = 0; // [-50:1:30]
right_ankle_roll = 0; // [-15:1:15]
right_shoulder_pitch = 0; // [-177:1:153]
right_shoulder_roll = 0; // [-129:1:91]
right_shoulder_yaw = 0; // [-150:1:150]
right_elbow = 0; // [-60:1:120]
/* [Display] */
show_lids = true;
/* [Hidden] */
c_white = [0.92,0.92,0.90]; c_lid = [0.85,0.86,0.88]; c_dark = [0.30,0.31,0.34];
 color(c_white) import("STL/assembly/pelvis_body.stl", convexity=10);
 if (show_lids) color(c_lid) import("STL/assembly/pelvis_lid1.stl", convexity=10);
 if (show_lids) color(c_lid) import("STL/assembly/pelvis_lid2.stl", convexity=10);
 if (show_lids) color(c_lid) import("STL/assembly/pelvis_lid3.stl", convexity=10);
translate([0.0000,22.5582,241.6050]) rotate(a=left_hip_pitch, v=[0.0000,1.0000,0.0000]) translate([-0.0000,-22.5582,-241.6050]) {
   color(c_white) import("STL/assembly/left_hip_pitch_link_body.stl", convexity=10);
   if (show_lids) color(c_lid) import("STL/assembly/left_hip_pitch_link_lid1.stl", convexity=10);
  translate([0.0000,40.7582,230.9423]) rotate(a=left_hip_roll, v=[0.9847,0.0000,0.1740]) translate([-0.0000,-40.7582,-230.9423]) {
     color(c_white) import("STL/assembly/left_hip_roll_link_body.stl", convexity=10);
    translate([16.1762,40.7582,189.6856]) rotate(a=left_hip_yaw, v=[-0.1740,0.0000,0.9847]) translate([-16.1762,-40.7582,-189.6856]) {
       color(c_white) import("STL/assembly/left_hip_yaw_link_body.stl", convexity=10);
       if (show_lids) color(c_lid) import("STL/assembly/left_hip_yaw_link_lid1.stl", convexity=10);
       if (show_lids) color(c_lid) import("STL/assembly/left_hip_yaw_link_lid2.stl", convexity=10);
      translate([-0.0008,41.5103,123.7965]) rotate(a=left_knee, v=[0.0000,1.0000,0.0000]) translate([0.0008,-41.5103,-123.7965]) {
         color(c_white) import("STL/assembly/left_knee_link_body.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/left_knee_link_lid1.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/left_knee_link_lid2.stl", convexity=10);
        translate([-0.0008,41.4773,18.7930]) rotate(a=left_ankle_pitch, v=[0.0000,1.0000,0.0000]) translate([0.0008,-41.4773,-18.7930]) {
           color(c_dark) import("STL/assembly/left_ankle_pitch_link_body.stl", convexity=10);
          translate([-0.0008,41.4773,12.6477]) rotate(a=left_ankle_roll, v=[1.0000,0.0000,-0.0000]) translate([0.0008,-41.4773,-12.6477]) {
             color(c_white) import("STL/assembly/left_ankle_roll_link_body.stl", convexity=10);
          }
        }
      }
    }
  }
}
translate([0.0000,-22.5582,241.6050]) rotate(a=right_hip_pitch, v=[0.0000,1.0000,0.0000]) translate([-0.0000,22.5582,-241.6050]) {
   color(c_white) import("STL/assembly/right_hip_pitch_link_body.stl", convexity=10);
   if (show_lids) color(c_lid) import("STL/assembly/right_hip_pitch_link_lid1.stl", convexity=10);
  translate([0.0000,-40.7582,230.9423]) rotate(a=right_hip_roll, v=[0.9847,0.0000,0.1740]) translate([-0.0000,40.7582,-230.9423]) {
     color(c_white) import("STL/assembly/right_hip_roll_link_body.stl", convexity=10);
    translate([16.1762,-40.7582,189.6856]) rotate(a=right_hip_yaw, v=[-0.1740,0.0000,0.9847]) translate([-16.1762,40.7582,-189.6856]) {
       color(c_white) import("STL/assembly/right_hip_yaw_link_body.stl", convexity=10);
       if (show_lids) color(c_lid) import("STL/assembly/right_hip_yaw_link_lid1.stl", convexity=10);
       if (show_lids) color(c_lid) import("STL/assembly/right_hip_yaw_link_lid2.stl", convexity=10);
      translate([-0.0008,-41.5103,123.7965]) rotate(a=right_knee, v=[0.0000,1.0000,0.0000]) translate([0.0008,41.5103,-123.7965]) {
         color(c_white) import("STL/assembly/right_knee_link_body.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/right_knee_link_lid1.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/right_knee_link_lid2.stl", convexity=10);
        translate([-0.0008,-41.4773,18.7930]) rotate(a=right_ankle_pitch, v=[0.0000,1.0000,0.0000]) translate([0.0008,41.4773,-18.7930]) {
           color(c_dark) import("STL/assembly/right_ankle_pitch_link_body.stl", convexity=10);
          translate([-0.0008,-41.4773,12.6477]) rotate(a=right_ankle_roll, v=[1.0000,0.0000,-0.0000]) translate([0.0008,41.4773,-12.6477]) {
             color(c_white) import("STL/assembly/right_ankle_roll_link_body.stl", convexity=10);
          }
        }
      }
    }
  }
}
translate([-1.3872,0.0000,292.9500]) rotate(a=waist_yaw, v=[0.0000,0.0000,1.0000]) translate([1.3872,-0.0000,-292.9500]) {
   color(c_white) import("STL/assembly/torso_body.stl", convexity=10);
  translate([-0.0025,35.0770,379.6730]) rotate(a=left_shoulder_pitch, v=[0.0002,0.9612,0.2757]) translate([0.0025,-35.0770,-379.6730]) {
     color(c_white) import("STL/assembly/left_shoulder_pitch_link_body.stl", convexity=10);
     if (show_lids) color(c_lid) import("STL/assembly/left_shoulder_pitch_link_lid1.stl", convexity=10);
     if (show_lids) color(c_lid) import("STL/assembly/left_shoulder_pitch_link_lid2.stl", convexity=10);
    translate([0.0001,49.1962,378.6865]) rotate(a=left_shoulder_roll, v=[1.0000,-0.0002,-0.0001]) translate([-0.0001,-49.1962,-378.6865]) {
       color(c_white) import("STL/assembly/left_shoulder_roll_link_body.stl", convexity=10);
      translate([-0.0014,51.3823,342.5666]) rotate(a=left_shoulder_yaw, v=[0.0001,-0.0001,1.0000]) translate([0.0014,-51.3823,-342.5666]) {
         color(c_white) import("STL/assembly/left_shoulder_yaw_link_body.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/left_shoulder_yaw_link_lid1.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/left_shoulder_yaw_link_lid2.stl", convexity=10);
        translate([5.5211,51.3830,314.3850]) rotate(a=left_elbow, v=[0.0002,1.0000,0.0001]) translate([-5.5211,-51.3830,-314.3850]) {
           color(c_white) import("STL/assembly/left_elbow_body.stl", convexity=10);
        }
      }
    }
  }
  translate([-0.0025,-35.0735,379.6730]) rotate(a=right_shoulder_pitch, v=[-0.0002,0.9612,-0.2757]) translate([0.0025,35.0735,-379.6730]) {
     color(c_white) import("STL/assembly/right_shoulder_pitch_link_body.stl", convexity=10);
     if (show_lids) color(c_lid) import("STL/assembly/right_shoulder_pitch_link_lid1.stl", convexity=10);
     if (show_lids) color(c_lid) import("STL/assembly/right_shoulder_pitch_link_lid2.stl", convexity=10);
    translate([0.0001,-49.1927,378.6865]) rotate(a=right_shoulder_roll, v=[1.0000,0.0002,-0.0001]) translate([-0.0001,49.1927,-378.6865]) {
       color(c_white) import("STL/assembly/right_shoulder_roll_link_body.stl", convexity=10);
      translate([-0.0014,-51.3788,342.5666]) rotate(a=right_shoulder_yaw, v=[0.0001,0.0001,1.0000]) translate([0.0014,51.3788,-342.5666]) {
         color(c_white) import("STL/assembly/right_shoulder_yaw_link_body.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/right_shoulder_yaw_link_lid1.stl", convexity=10);
         if (show_lids) color(c_lid) import("STL/assembly/right_shoulder_yaw_link_lid2.stl", convexity=10);
        translate([5.5211,-51.3795,314.3850]) rotate(a=right_elbow, v=[-0.0002,1.0000,-0.0001]) translate([-5.5211,51.3795,-314.3850]) {
           color(c_white) import("STL/assembly/right_elbow_body.stl", convexity=10);
        }
      }
    }
  }
   color(c_dark) import("STL/assembly/head_body.stl", convexity=10);
  color([0.95,0.75,0.2]) import("STL/assembly/head_brow_left.stl");
  color([0.95,0.75,0.2]) import("STL/assembly/head_brow_right.stl");
}
