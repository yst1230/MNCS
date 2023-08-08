clc;clear
name_lists={'aa','al','av','aw','ay'};
for i=1:5
    load(['data_set_IVa_' name_lists{i} '.mat'])
    load(['true_labels_' name_lists{i} '.mat'])
    event_imagery=mrk.pos;
    for j=1:length(true_y)
        EEGSignals(j,:,:)=0.1*double(cnt(event_imagery(j)+1:event_imagery(j)+3.5*nfo.fs,:));
    end
    save(['eegdata_' name_lists{i}],'EEGSignals')
end