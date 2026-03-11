fs = 44100; recObj = audiorecorder(fs,16,1);
disp('Say the letter'); recordblocking(recObj, 2);
y = getaudiodata(recObj); audiowrite('a.wav', y, fs);