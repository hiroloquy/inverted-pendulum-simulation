reset

#=================== Parameter ====================
# Wheel
wid_w = 0.1 	# Width [m]
dens_w = 1000	# Density [kg/m3]
radius_w = 0.5  # Radius [m]
mass_w = (pi * radius_w ** 2) * wid_w * dens_w	# Mass [kg]

# Inverted pendulum
wid_p = 0.1		# Width [m]
dens_p = 1000	# Density [kg/m3]
len_p= 2.5		# Length [m]
mass_p = wid_p * 2 * len_p* dens_p # Mass [kg]

dist_cen_grav = 1.0	# Distance between the centers of rotation and gravity of the inverted pendulum [m]
dpl = dist_cen_grav + len_p/2
lmd = len_p/2 - dist_cen_grav

# Gains
# Successful example
kp1 = 10000.0	# Proportional gain of the inverted pendulum [Nm/rad] (theta)
kp2 = 300.0   # Proportional gain of the wheel [Nm/rad] (phi)
kd1 = 2500.0   # Differential gain of the inverted pendulum [Nms/rad] (dtheta)
kd2 = 500.0   # Differential gain of the wheel [Nms/rad] (dphi)

# Failure example (vibration)
# kp1 = 12000.0	# Proportional gain of the inverted pendulum [Nm/rad] (theta)
# kp2 = 300.0   # Proportional gain of the wheel [Nm/rad] (phi)
# kd1 = 2500.0   # Differential gain of the inverted pendulum [Nms/rad] (dtheta)
# kd2 = 600.0   # Differential gain of the wheel [Nms/rad] (dphi)
# kd2 = 700.0   # Differential gain of the wheel [Nms/rad] (dphi)

# Simulation Conditions
dt = 0.002	        # Time step [s]
dh = dt/6	        # Coefficient of Runge Kutta 4th [s]
g = 9.81	        # Gravitational acceleration [m/s2]
TIME_MAX = 10.0	    # Maximum of time[s] 
STEP_MAX = TIME_MAX/dt
NUM_SKIP = 5        # Number of the skipping frame
DELAY_TIME = 0.0001 # [s] (This parameter is valid only in qt terminal.)
plotRange = 1.5

# Initial values
t  = 0.0			    # [s]
# Pendulum
theta = pi/180*10 	    # [rad]
dtheta = 0              # [rad/s]
# Wheel
phi = 0.0			    # [rad]
dphi = 0.0			    # [rad/s]

# Select terminal type
qtMode = 0     # ==1: qt (simulator) / !=1: png (output images for making video)
print sprintf("[MODE] %s", (qtMode==1 ? 'Simulate in Qt window' :'Output PNG images'))

#=================== Function ====================
# Torque generated from the wheel's motor
u(theta, dtheta, phi, dphi) = kp1*theta + kd1*dtheta + kp2*phi + kd2*dphi

coef_mag = mass_p * dist_cen_grav * g
coef_mar = mass_p * dist_cen_grav * radius_w

# Matrix A
A11 = mass_p * (dist_cen_grav ** 2 + (len_p ** 2) / 12.)
A12(theta)= coef_mar * cos(theta)
A22 = (mass_p + 3./2. * mass_w) * radius_w ** 2
detA(theta) = A11*A22 - A12(theta)**2

B1(theta, dtheta, phi, dphi) = - u(theta, dtheta, phi, dphi) + coef_mag*sin(theta) 
B2(theta, dtheta, phi, dphi) = u(theta, dtheta, phi, dphi) + coef_mar*sin(theta)*(dtheta)**2

# ODE
# d??/dt
f1(theta, dtheta, phi, dphi) = dtheta
# d2??/dt
f2(theta, dtheta, phi, dphi) = \
  (A22*B1(theta, dtheta, phi, dphi) - A12(theta)*B2(theta, dtheta, phi, dphi)) / detA(theta)
# d??/dt
f3(theta, dtheta, phi, dphi) = dphi
# d2??/dt
f4(theta, dtheta, phi, dphi) = \
  (A11*B2(theta, dtheta, phi, dphi) - A12(theta)*B1(theta, dtheta, phi, dphi)) / detA(theta)

# Runge-Kutta 4th (Define rk_i(theta, dtheta, t))
do for[i=1:4]{
    rki = "rk"
    fi  = "f".sprintf("%d", i)
    rki = rki.sprintf("%d(theta, dtheta, phi, dphi) = (\
        k1 = %s(theta, dtheta, phi, dphi),\
        k2 = %s(theta + dt*k1/2., dtheta + dt*k1/2., phi + dt*k1/2., dphi + dt*k1/2.),\
        k3 = %s(theta + dt*k2/2., dtheta + dt*k2/2., phi + dt*k2/2., dphi + dt*k2/2.),\
        k4 = %s(theta + dt*k3, dtheta + dt*k3, phi + dt*k3, dphi + dt*k3),\
        dh * (k1 + 2*k2 + 2*k3 + k4))", i, fi, fi, fi, fi)
    eval rki
}

# Show the value of time t
showTime(t) = sprintf("{/:Italic t} = %2.2f s", t) 

#=================== Settimg Terminal ====================
if(qtMode==1){
    set term qt size 1280, 720 font 'Times'
} else {
    set term pngcairo size 1280, 720 font 'Times'
    folderName = 'png_multiplot'
    system sprintf('mkdir %s', folderName)
}

#=================== Calculation ====================
outputData = "data.txt"
print sprintf('Start outputting %s ...', outputData)
set print outputData

# Write terms and initial values into txt file
print sprintf("# %s", outputData)
print sprintf("# dt=%.5f, g=%.5f", dt, g)
print "# t / ?? / ??' / ?? / ??'"
print sprintf("%.3f %.3f %.3f %.3f %.3f", t, theta, dtheta, phi, dphi)

do for [i=1:STEP_MAX:1] {
    t  = t  + dt

    incTheta = rk1(theta, dtheta, phi, dphi)
    incDtheta = rk2(theta, dtheta, phi, dphi)
    incPhi = rk3(theta, dtheta, phi, dphi)
    incDphi = rk4(theta, dtheta, phi, dphi)

    theta = theta + incTheta
    dtheta = dtheta + incDtheta
    phi = phi + incPhi
    dphi = dphi + incDphi
	
	print sprintf("%.3f %.3f %.3f %.3f %.3f", t, theta, dtheta, phi, dphi)
}

unset print 
print sprintf('Finish!')

 #=================== Plot ====================
if(qtMode == 1) {
    print "Start simulation"
} else {
    print sprintf('Start outputting %d images ...', STEP_MAX/NUM_SKIP+1)
}

do for [i=0:int(STEP_MAX/NUM_SKIP):1] {
    if(qtMode != 1) { set output sprintf("%s/img_%04d.png", folderName, i) }

    # Get the value of time from either of txt files
    get_data_num = i*NUM_SKIP
    set yrange [*:*]    # This command enables to remove restrictions on the range of the stats command.
    stats outputData using 1 every ::get_data_num::get_data_num nooutput; t = STATS_max
    stats outputData using 2 every ::get_data_num::get_data_num nooutput; theta = STATS_max
    stats outputData using 4 every ::get_data_num::get_data_num nooutput; phi = STATS_max

	set multiplot
	# Left side: Segway animation
        unset key
        set grid
        set xlabel '{/:Italic x} [m]' font ', 20'
        set ylabel '{/:Italic y} [m]' font ', 20' offset screen -0.015, 0
        set tics 0.5 font ', 18'
        set size ratio -1 0.9
        set origin -0.4, 0.03

        xg =radius_w * phi	# X coordinate of the wheel center
        set object 1 circle at xg, radius_w fc rgb 'gray50' size radius_w fs solid behind # Wheel
        set object 2 circle at xg, radius_w fc rgb 'black' size 0.03 fs solid front	# Point representing the center of the wheel
        set arrow 1 nohead from xg, radius_w to xg+radius_w*sin(phi), radius_w*(1+cos(phi)) lc -1 lw 3 front # Radius of the wheel
        set arrow 2 nohead from xg-lmd*sin(theta), radius_w-lmd*cos(theta) to xg+dpl*sin(theta), radius_w+dpl*cos(theta) lc 6 lw 12 back # Pendulum
        set label 1 left showTime(t) at graph 0.4, 1.05 font ', 20'
        plot [xg-plotRange:xg+plotRange][-0.5:3.5] 0 lc -1 lw 2 # Ground
        
        # Unset objects, arrows, and labels in order not to display them in the other plot area
        unset object
        unset arrow
        unset label

	# Right side: Graph
        set key right top font ', 20'
        set grid
        set xl '{/:Italic t} [s]' font ', 20'
        set yl '{/:Italic ??, ??} [rad]' font ', 20' offset screen -0.015, 0
        # set yl '{/:Italic ??, ??} [rad]  / {/:Italic d??, d??} [rad/s] ' font ', 20' offset screen -0.03, 0
        set xr[0:TIME_MAX]
        set yr[-3:3]
        set tics 1 font ', 18'
        set origin 0.3, 0.03
        set size ratio 1.1

        plot outputData using 1:2 every ::0::get_data_num with lines lw 3 lc 4 t "{/:Italic ??}  ",\
            outputData using 1:4 every ::0::get_data_num with lines lw 3 lc 7 t "{/:Italic ??}  ",\
            0 lw 1 lc -1 notitle

        # plot outputData using 1:2 every ::0::get_data_num with lines lw 3 lc 4 t "{/:Italic ??}  ",\
        #      outputData using 1:3 every ::0::get_data_num with lines lw 3 lc 6 t "{/:Italic d??} ",\
        #      outputData using 1:4 every ::0::get_data_num with lines lw 3 lc 7 t "{/:Italic ??}  ",\
        #      outputData using 1:5 every ::0::get_data_num with lines lw 3 lc 1 t "{/:Italic d??} ",\
        #      0 lw 1 lc -1 notitle
	unset multiplot

    if(qtMode == 1) {   
        pause ((i == 0 || i == int(STEP_MAX/NUM_SKIP)) ? 1 : DELAY_TIME)   # Adjust the drawing speed [s]
    } else {
        set out     # terminal pngcairo
    }
}

set out
print sprintf('Finish!')
