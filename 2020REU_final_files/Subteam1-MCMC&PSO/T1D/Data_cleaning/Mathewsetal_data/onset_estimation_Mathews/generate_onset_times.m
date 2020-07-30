%This code takes the data extacted from Mathews et al (2010)
%using grabit on Figure 2(E) and fits lognormal distributions to the onset
%data. We next extract data from the fitted distributions and uses the
%extracted onset times to extract n=109 for acute and n=278 for
%progressive.


%B. Shtylla (07/31/2019)
% close all 
% clear all
%============================
% ACUTE MICE ================
%============================

% open the csv file for acute mice
acute_df = textscan( regexprep( fileread('acute-age-onset.csv'), '\$', '0' ), '%f%f', 'delimiter', ',', 'HeaderLines', 1);
% retrieve only the ages from the second row in days(multiply by 7)
acute_agesweeks = acute_df{2};
% get normal distribution from data
acute_distweeks = fitdist(acute_agesweeks, 'lognormal');

%Sample the fitted distribution to obtain the experimental numbers
sample_acuteweeks = random(acute_distweeks, 109,1);

acute_agesdays = 7*acute_df{2};

% get normal distribution from data
acute_distdays = fitdist(acute_agesdays, 'lognormal');

%Sample the fitted distribution to obtain the experimental numbers
sample_acutedays = random(acute_distdays, 109,1);

%Now generate a shifted acute data set
modified_acute=importdata('converted_raw_acute.csv');
modified_acute=modified_acute.data;
resample_acute=random(acute_distdays, length(modified_acute(:,1)),1);


for i=1:length(modified_acute)-1
shifted_acute(:,1+(i-1)*2)=modified_acute(1,:)+round(resample_acute(i));
shifted_acute(:,i*2)=modified_acute(i+1,:);

unshifted_acute(:,1+(i-1)*2)=modified_acute(1,:);
unshifted_acute(:,i*2)=modified_acute(i+1,:);

end


%============================
% PROGRESSIVE MICE ==========
%============================

% open the csv file for progressive mice
progressive_df = textscan( regexprep( fileread('progressive-age-onset.csv'), '\$', '0' ), '%f%f', 'delimiter', ',', 'HeaderLines', 1);
% retrieve only the ages from the second row

progressive_agesweeks = progressive_df{2};
% get normal distribution from data
progressive_distweeks = fitdist(progressive_agesweeks, 'lognormal');
%Sample the fitted distribution to obtain the experimental numbers
sample_progressiveweeks = random(progressive_distweeks, 278,1);

progressive_agesdays = 7*progressive_df{2};
% get normal distribution from data
progressive_distdays = fitdist(progressive_agesdays, 'lognormal');
%Sample the fitted distribution to obtain the experimental numbers
sample_progressivedays = random(progressive_distdays, 278,1);

%Now generate a shifted acute data set
modified_progressive=importdata('converted_raw_progressive.csv');
modified_progressive=modified_progressive.data;
resample_progressive=random(progressive_distdays, length(modified_progressive(:,1)),1);


for i=1:length(modified_progressive)-1
    
shifted_progressive(:,1+(i-1)*2)=modified_progressive(1,:)+round(resample_progressive(i));
shifted_progressive(:,i*2)=modified_progressive(i+1,:);

unshifted_progressive(:,1+(i-1)*2)=modified_progressive(1,:);
unshifted_progressive(:,i*2)=modified_progressive(i+1,:);
end

% 
% figure;
% subplot(3,2,1)
% % generate a histogram from age onset column, 15 bins, lognormal fit
% histfit(acute_agesweeks, 30, 'lognormal');
% % limit x-axis scale
% set(gca, 'Xlim', [8, 30]);
% % label x and y axises
% ylabel('Frequencies');
% xlabel('Day Of T1D Onset');
% title('Data Fit: Acute Age Onset', 'fontsize', 16);
% subplot(3,2,2)
% % generate a histogram from age onset column, 30 bins, lognormal fit
% progressive_his = histfit(progressive_agesweeks, 30, 'lognormal');
% % limit x-axis scale
% set(gca, 'Xlim', [9, 30]);
% % label x and y axises
% ylabel('Frequencies');
% xlabel('Day Of T1D Onset');
% title('Data Fit: Progressive Age Onset', 'fontsize', 16);
% subplot(3,2,3)
% % generate a histogram from age onset column, 30 bins, lognormal fit
% h2=histogram(sample_acuteweeks,100);
% hold on
% h3=histfit(sample_acuteweeks,100,'lognormal');
% hold off
% % limit x-axis scale
% set(gca, 'Xlim', [9, 30]);
% h2.FaceColor=[0 0.5 0.5];
% % label x and y axises
% ylabel('Frequencies');
% xlabel('Day Of T1D Onset');
% title('Our Sampling: Acute Age Onset', 'fontsize', 16);
% 
% subplot(3,2,4)
% % generate a histogram from age onset column, 30 bins, lognormal fit
% h2=histogram(sample_progressiveweeks,100);
% hold on
% h3=histfit(sample_progressiveweeks,100,'lognormal');
% hold off
% % limit x-axis scale
% set(gca, 'Xlim', [9, 30]);
% h2.FaceColor=[0 0.5 0.5];
% % label x and y axises
% ylabel('Frequencies');
% xlabel('Day Of T1D Onset');
% title('Our Sampling: Progressive Age Onset', 'fontsize', 16);
% 
% subplot(3,2,5)
% for i=1:length(modified_acute)-1
% plot(shifted_acute(:,1+(i-1)*2),shifted_acute(:,i*2))
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Sampled Glucose Levels');
% title('Our Glucose: Acute Onset', 'fontsize', 16);
% 
% subplot(3,2,6)
% for i=1:length(modified_progressive)-1
% plot(shifted_progressive(:,1+(i-1)*2),shifted_progressive(:,i*2))
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Sampled Glucose Levels');
% title('Our Glucose: Progressive Onset', 'fontsize', 16);
% 
% 
% figure
% 
% subplot(1,2,1)
% for i=1:length(modified_acute)-1
% plot(unshifted_acute(:,1+(i-1)*2),shifted_acute(:,i*2),'k')
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Original Glucose Levels');
% title('Acute Onset', 'fontsize', 16);
% 
% subplot(1,2,2)
% for i=1:length(modified_progressive)-1
% plot(unshifted_progressive(:,1+(i-1)*2),shifted_progressive(:,i*2),'k')
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Original Glucose Levels');
% title('Progressive Onset', 'fontsize', 16);
% 
% figure
% 
% subplot(1,2,1)
% for i=1:length(modified_acute)-1
% plot(shifted_acute(:,1+(i-1)*2),shifted_acute(:,i*2))
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Our Sampled Glucose Levels');
% title('Acute Onset', 'fontsize', 16);
% 
% subplot(1,2,2)
% for i=1:length(modified_progressive)-1
% plot(shifted_progressive(:,1+(i-1)*2),shifted_progressive(:,i*2))
% hold all
% end
% hold off
% xlabel('Time (days)');
% ylabel('Our Sampled Glucose Levels');
% title('Progressive Onset', 'fontsize', 16);

