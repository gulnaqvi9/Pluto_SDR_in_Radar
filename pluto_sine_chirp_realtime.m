clear;
clc;

initradar();
load('radar_init_data_sine_chirp.mat')
[txW]=createwaveforms(fs,nsinesamp,chirpsamples,nblank,win_tx_pulse);

readfrompluto=1;
if(readfrompluto)
    [tx,rx]=initipluto(1e9,fs);
    [rawdata]=fastdatacollect(tx,rx,txW);
else
    load('rawdata.mat')
    rawdata=rdsyncret;
end
processdata(rawdata(:,1:chirpsamples+nblank));


%% process data
function processdata(rdsync)
    load 'radar_init_data_sine_chirp.mat' fstart
    load 'radar_init_data_sine_chirp.mat' fend
    load 'radar_init_data_sine_chirp.mat' c
    load 'radar_init_data_sine_chirp.mat' fs
    load 'radar_init_data_sine_chirp.mat' chirpsamples
    [nfreq,framelength]=size(rdsync);
    
    [txWchirp]=createchirp(fs,chirpsamples,2);
    
    xx=xcorr2(rdsync,txWchirp.');
    rdcorr=xx(:,1:framelength);
    figure
    [pki,lci]=findpeaks(abs(rdcorr(1,:)));
    %size(pki)
    plot(20*log10(abs(rdcorr(1,:))));
    hold on
    plot(lci,20*log10(abs(pki)),'x');
    getframe();

    Tmajstep=1/fs;
    %dmajstep=c*Tmajstep/2;
    time=0:Tmajstep:(framelength-1)*Tmajstep;
    dmaj=time*c/2;
    %ntime=length(time);


    nfft=2^(nextpow2(8*nfreq));
    xst=nfft/nfreq;
    range= c*(0:nfft-1)/(fend-fstart)/2/xst;

    for jj=1:nfreq
        %[pki,lci]=findpeaks(abs(rdcorrsave(jj,:)));
        corrpeaks(jj,1:length(lci))=rdcorr(jj,lci);
        locs(jj,1:length(lci))=lci;
        dlocs(jj,1:length(lci))=dmaj(lci);
    end
    hwinf=hanning(nfreq);
    corrpeaksf=corrpeaks.*hwinf;
    signalfft=ifft(corrpeaksf,nfft);
    signalfft=signalfft/max(max(abs(signalfft)));
    [nf,npeaks]=size(corrpeaks);
    figure(3)
    for jj=1:npeaks
        plot(range,20*log10(abs(signalfft(:,jj))),'LineWidth',1);
        xlabel('distance (m)');
        ylabel('Amplitude');
        ylim([-60 10]);
        xlim([0 range(end)]);
        grid on;
        hold on
        getframe();
    end
    xlabel('distance (m)');
    ylabel('Amplitude');

    % gate the first range

   % td=signalfft(1:22,7);
   % figure
   % plot(unwrap(angle(td)))
    
end
%% initialize radar parameters
function initradar()
    c=3e8;

    fs = 60e6;  
    nsinesamp=100;
    
    chirpsamples=100;
    win_tx_pulse=2;
    
    nblank=300;

    % Frequency range for Fourier Transformation
    fstart=1.2e9;
    fend=1.4e9;
    fstep=2e6;

    save('radar_init_data_sine_chirp')
end



%% initialize pluto
function [tx,rx]=initipluto(freq,fs)
    tx = sdrtx('Pluto');
    tx.CenterFrequency = freq;
    tx.BasebandSampleRate = fs;
    tx.Gain = 0;
    
    rx=sdrrx('Pluto');
    rx.CenterFrequency = freq;
    rx.BasebandSampleRate = fs;
    % rx.GainSource='Manual';
    % rx.Gain=70;
end

%% Fast data collect

function [rdsyncret]=fastdatacollect(tx,rx,txW)
    load 'radar_init_data_sine_chirp.mat' fstart
    load 'radar_init_data_sine_chirp.mat' fend
    load 'radar_init_data_sine_chirp.mat' fstep
    freq=fstart:fstep:fend;
    framelength=length(txW);
    nfreq=length(freq);
   % nfreq=5;
    rdsyncret=zeros(nfreq,framelength);
    tic

  
    
    for ii=1:nfreq
        fprintf('Processing freq: %8.5f \n',freq(ii)/1e9);

        tx.CenterFrequency = freq(ii);
        transmitRepeat(tx,txW);
        rx.CenterFrequency = freq(ii);
      
        rd=capture(rx,.01,'Seconds');

 

        %release(rx);
       % release(tx);



        rd=cast(rd,"double");
    
        [rdmean]=meanoverframe(rd,framelength);
    
%         figure(ii+1)
%         plot(real(rdmean),'r-','LineWidth',1)
%         hold on
%         plot(imag(rdmean),'b-','LineWidth',1)
        %plot(abs(rdmean))
    
        [rdsync]=syncwaveforms(rdmean);
    
        phs=angle(rdsync(11));
        %rdsync=rdsync/rdsync(11);
        rdsync=rdsync*exp(-1j*phs);
        
        rdsyncret(ii,:)=rdsync.';
        
%         figure(1)
%         plot(real(rdsync),'r-','LineWidth',1)
%         hold on
%         plot(imag(rdsync),'b-','LineWidth',1)
%         plot(abs(rdsync))
%         xlim([1 2*framelength]);
%     
%         getframe();

    end
    toc
    release(tx)
    release(rx)
    save('rawdata')
end
%% collect data from pluto
function [rdsyncret]=datacollect(tx,rx,txW)
    load 'radar_init_data_sine_chirp.mat' fstart
    load 'radar_init_data_sine_chirp.mat' fend
    load 'radar_init_data_sine_chirp.mat' fstep
    freq=fstart:fstep:fend;
    framelength=length(txW);
    nfreq=length(freq);
    nfreq=5;
    rdsyncret=zeros(nfreq,framelength);
    tic
    for ii=1:nfreq
        fprintf('Processing freq: %8.5f \n',freq(ii)/1e9);

        tx.CenterFrequency = freq(ii);
        transmitRepeat(tx,txW);
        rx.CenterFrequency = freq(ii);
      
        rd=capture(rx,.05,'Seconds');
       
 

        release(rx);
        release(tx);



        rd=cast(rd,"double");
    
        [rdmean]=meanoverframe(rd,framelength);
    
%         figure(ii+1)
%         plot(real(rdmean),'r-','LineWidth',1)
%         hold on
%         plot(imag(rdmean),'b-','LineWidth',1)
        %plot(abs(rdmean))
    
        [rdsync]=syncwaveforms(rdmean);
    
        phs=angle(rdsync(11));
        %rdsync=rdsync/rdsync(11);
        rdsync=rdsync*exp(-1j*phs);
        
        rdsyncret(ii,:)=rdsync.';
        
%         figure(1)
%         plot(real(rdsync),'r-','LineWidth',1)
%         hold on
%         plot(imag(rdsync),'b-','LineWidth',1)
%         plot(abs(rdsync))
%         xlim([1 2*framelength]);
%     
%         getframe();

    end
    toc
    save('rawdata')
end


%% average the waveforms
function [rdmean]=meanoverframe(rd,framelength)
    nsamples=fix(length(rd)/framelength);
    rdsum=zeros(framelength,1);
%     figure(7)
%     hold on;
    for jj=1:nsamples
        rdsum=rdsum+rd((jj-1)*framelength+1:(jj)*framelength);
%         plot(abs(rd((jj-1)*framelength+1:(jj)*framelength)));
%         getframe;
    end
    rdmean=rdsum/nsamples;
   %plot(abs(rdmean),'b');
end

%% Function to sync the waveforms
function [rdsync]=syncwaveforms(rdmean)
load 'radar_init_data_sine_chirp.mat' fs
 load 'radar_init_data_sine_chirp.mat' chirpsamples

    %normratio=max(abs(txW))/max(abs(rdmean));
    %txW1=txW/normratio;
    %refframe=circshift(abs(txW1),10);
    %rdsave(ii,:)=circshift(rdsave(ii,:),fix(chirpsamples/2));
    %e1=abs(rdmean);
    %plot(e1,'r');
    %hold on
    %plot(refframe,'r');
   % try
        %[~,~,D] = alignsignals(refframe,e1,Method="risetime",StateLevels=[100 600]);

        [txWchirp]=createchirp(fs,chirpsamples,1);
        [~,~,D] = alignsignals(txWchirp,rdmean,Method="xcorr");
   % catch
   %     fprintf('error in alignsignals \n')
        
   % end
%     if (D <= 0)
%         fprintf('%d \n', D)
%         return
%     end
    rdsync=circshift(rdmean,-D);
end
%% function to create sinusoidal signal
function [txWsine]=createsine(fs,nsinesamp)
    sw = dsp.SineWave;
    sw.Amplitude = 0.1;
    sw.Frequency = 1e6;
    sw.PhaseOffset=0;
    sw.ComplexOutput = true;
    sw.SampleRate = fs;
    sw.SamplesPerFrame = nsinesamp;
    txWsine = sw();  
end


%% function to create chirp signal
function [txWchirp]=createchirp(fs,chirpsamples,win_tx_pulse)
    tchirp=0:1/fs:(chirpsamples-1)/fs;
    nchirp=length(tchirp);
    twin=tukeywin(nchirp,0.2);
    hwin=hanning(nchirp);
    uwin=ones(nchirp,1);
    BW=10e6;
    fc=5.01e6;
    txW1=chirp(tchirp,(fc-BW/2),chirpsamples/fs,(fc+BW/2),'complex');
    txWchirp=txW1.';
    if (win_tx_pulse == 0)
        txWchirp=uwin.*txWchirp;
    elseif (win_tx_pulse ==1)
        txWchirp=twin.*txWchirp;
    elseif (win_tx_pulse==2)
        txWchirp=hwin.*txWchirp;
    else
        fprintf('win_tx_pulse should be 0 through 2 only and cannot be %d \n',win_tx_pulse);
        return;
    end

end

%% create waveforms
function [txW]=createwaveforms(fs,nsinesamp,chirpsamples,nblank,win_tx_pulse)
    [txWsine]=createsine(fs,nsinesamp);
    [txWchirp]=createchirp(fs,chirpsamples,win_tx_pulse);
    blanksig=zeros(nblank,1);
    framelength=nsinesamp+nblank+chirpsamples+nblank;
    txW=[txWsine; blanksig; txWchirp; zeros(nblank,1)];
end