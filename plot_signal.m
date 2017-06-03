function []  = plot_signal( data, subject,trial_start, trial_end )
    figure
    for k=1:25
        signal_type='EEG';
        color='blue';
        if k>22
            signal_type = 'EOG';
            color='red';
        end
        subplot(5,5,k);
        plot([trial_start:trial_end],data{1,subject}.X(trial_start:trial_end,k:k),color);
        xlim([trial_start trial_end]);
        title(sprintf('%s channel %s',signal_type,int2str(k)));
    end
end

