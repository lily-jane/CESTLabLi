function satpulse = getsatpulse(curve, angle, t, duration, TR)
num = 360*sum(abs(curve))*duration/(256*angle*2*pi); 
satpulse = curve./num;
satpulse = satpulse(t);
end