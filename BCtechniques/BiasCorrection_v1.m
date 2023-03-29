%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bias correction for simulation data via observation data
%%% Markel Penalba, Mondragon Unibertsitatea
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc

%% Definition of the case study and data loading
%%%
disp('Loading data... \n') ;
%%%
% Observation data (Puertos del Estado - Spanish Oceanography agency: https://www.puertos.es/es-es/oceanografia/Paginas/portus.aspx)
Dir_Obs='C:\Users\mpenalba\Dropbox sMGEP)\Ikerketa\0_MREresource\0_MetoceanData\MeasuredData';
FolderName_Obs='PuestosEstado' ;
FileName_Obs='CaboBegur' ;        % {'BilbaoVizcaya','CaboSilleiro','GolfoCadiz','CaboBegur','GranCanaria'} ;
TimePeriod_Obs=[2000 2019] ;
    FullFileName_Obs=[Dir_Obs,'\',FolderName_Obs,'-',FileName_Obs,'_',...
        num2str(TimePeriod_Obs(1)),'-',num2str(TimePeriod_Obs(end)),'_clean.csv'] ;
    %%%
    MetoceanData0_Obs=readtable(FullFileName_Obs) ;
    %%%
% Simulation data (ECMWF-ERA5: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form)
Dir_Sim='C:\Users\mpenalba\Dropbox (MGEP)\Ikerketa\0_MREresource\0_MetoceanData\ReanalysisData\ERA5\csvFiles' ;
FolderName_Sim='CapCreus' ;     %{'GulfBiscay','CaboSilleiro','GulfCadiz','CapCreus','GranCanaria'} ;
FileName_Sim='MedFloat' ;     %{'SaitecGeroa','WindFloat2','MarDeAgata','MedFloat','Gofio'} ;
TimePeriod_Sim=[1995 2019] ;
    FullFileName_Sim=[Dir_Sim,'\',FolderName_Sim,'-',FileName_Sim,'_',...
        num2str(TimePeriod_Sim(1)),'-',num2str(TimePeriod_Sim(end)),'.csv'] ;
    %%%
    MetoceanData0_Sim=readtable(FullFileName_Sim) ;
    %%%

%% Data selection and preparation for bias correction
%%%
pause(1) ; clc
disp('Data preparation... \n') ;
%%%
% Selection of the time period for applying BC
SaveData=1 ;
Years0=[2000 2019] ;
while Years0(1)<min(year(table2array(MetoceanData0_Obs(:,1)))) || Years0(1)<min(year(table2array(MetoceanData0_Sim(:,1))))
    Years0(1)=Years0(1)+1 ;
end
while Years0(end)>max(year(table2array(MetoceanData0_Obs(:,1)))) || Years0(end)>max(year(table2array(MetoceanData0_Sim(:,1))))
    Years0(end)=Years0(end)-1 ;
end
ResultsFolderName=['Results/Figures/',FileName_Obs,'_',num2str(Years0(1)),'-',num2str(Years0(end))] ;
% Selection of the data
VarNames_var={'Hs [m]','Tp [s]','Uw [m/s]'} ;
if isa(table2array(MetoceanData0_Obs(1,1)),'datetime')
    Hs_column(1)=2 ; Tp_column(1)=3 ; Uw_column(1)=9 ;
else
    disp('Define the position of the requried variables! \n')
    return
end
if isa(table2array(MetoceanData0_Sim(1,1)),'datetime')
    Hs_column(2)=2 ; Tp_column(2)=3 ; Uw_column(2)=5 ;
else
    Hs_column(2)=5 ; Tp_column(2)=6 ; Uw_column(2)=8 ;
end
[MetoceanData_Obs,MetoceanData_Sim,InputSet_Datetime,InputSet_Obs,InputSet_Sim,h1,h2]=...
    DataPreparation(MetoceanData0_Obs,MetoceanData0_Sim,Years0,Hs_column,Tp_column,Uw_column) ;
    if ~exist(ResultsFolderName,'dir')
        mkdir(ResultsFolderName)
    end
    if SaveData
        set(h1,'units','normalized','outerposition',[0 0 1 1]) ;
        saveas(h1,[ResultsFolderName,'/TimeSeries'],'fig') ;
        saveas(h1,[ResultsFolderName,'/TimeSeries'],'svg') ;
        saveas(h2,[ResultsFolderName,'/PDFs'],'fig') ;
        saveas(h2,[ResultsFolderName,'/PDFs'],'svg') ;
    end

pause(1)
MetoceanData_Obs0=MetoceanData0_Obs(:,[1 Hs_column(1) Tp_column(1) Uw_column(1)]) ;
MetoceanData_Sim0=MetoceanData0_Sim(:,[1 Hs_column(2) Tp_column(2) Uw_column(2)]) ;

%% Bias Correction
clc
clearvars -except FullFileName_Obs FullFileName_Sim ResultsFolderName SaveData...
    MetoceanData_Obs0 MetoceanData_Sim0 Years0...
    MetoceanData_Obs MetoceanData_Sim VarNames_var...
    InputSet_Datetime InputSet_Obs InputSet_Sim
BCstats=cell(length(InputSet_Obs(1,:)),1) ;
polyParams_EQM2=cell(length(InputSet_Obs(1,:)),1) ;
Data_Sim_corrected=cell(length(InputSet_Obs(1,:)),1) ;
for di=1:length(InputSet_Obs(1,:))
    Datetime_var=InputSet_Datetime{1,di} ;
    Data_Obs=InputSet_Obs{1,di} ;
    Data_Sim=InputSet_Sim{1,di} ;
    qvec=[5:10:95 98 99 99.5 99.75 99.9 9.95 99.96 99.97 99.98 99.99] ;                 % Vector of the quantiles distribution
    % DELTA bias-correction technique
    BCd_Hs=Delta_BC(Data_Obs,Data_Sim) ;
    Data_Sim_Delta=Data_Sim+BCd_Hs ;

    % EQM bias-correction technique (#1 - full dataset)
    Prob_Data_Obs1=Prob_data(Data_Obs)' ; Prob_Data_Sim1=Prob_data(Data_Sim)' ;
    nPoly=3 ;
    [~,p_EQM1,Data_Sim_EQM1,bias_EQM1]=...
        f_CFDmatching(nPoly,Data_Obs,Prob_Data_Obs1,Data_Sim,Prob_Data_Sim1,0) ;

    % EQM bias-correction technique (#2 - Linearly spaced quantiles)
    qNum=50 ;                                       % Number of quantiles included in the BC
    nPoly=1 ;                                       % Polynomial order
    q_Obs(:,di)=quantile(Data_Obs,linspace(0.01,0.999,qNum))' ;          % qNum quantiles of observation data
    q_Sim=quantile(Data_Sim,qNum)' ;                % qNum quantiles of simulated data
    Data_Sim_EQM2=zeros(length(Data_Sim),1) ;
    [Data_Sim_EQM2_qi,pos_Ref,p_EQM2]=DataClass_quantiles_Obs(qNum,nPoly,...
        q_Obs(:,di),Data_Obs,q_Sim,Data_Sim,0) ;
% % %     [Data_Sim_EQM2_qi,pos_Ref,p_EQM2]=DataClass_quantiles_Sim(qNum,nPoly,...
% % %         q_Obs(:,di),Data_Obs,q_Sim,Data_Sim,0) ;
    posTOT=0 ; DataTOT=0 ;
    for qi=1:length(q_Obs(:,di))+1
        Data_Sim_EQM2(pos_Ref{qi,1})=Data_Sim_EQM2_qi{qi,1} ;
            posTOT=posTOT+length(pos_Ref{qi,1}) ;
            DataTOT=DataTOT+length(Data_Sim_EQM2_qi{qi,1}) ;
    end
        %%% Double checking if all points of the dataset have been incorporated
        if posTOT~=length(Data_Sim_EQM2) || DataTOT~=length(Data_Sim_EQM2)
            disp('Simulation data classiffication process is not correct! Lengths of the vector before and after are not consistent \n')
        end
% % %         figure;plot(Datetime_var,Data_Obs,'k',Datetime_var,Data_Sim,'b',Datetime_var,Data_Sim_EQM2,'r--',Datetime_var,Data_Sim_EQM2_2','g-.',Datetime_var,Data_Sim_EQM2_3','m:')

    % EGQM bias-correction technique
    qNum=49 ;                                       % Number of quantiles included in the BC
    nu=0 ; beta=1 ;
    x_qi=zeros(qNum,1) ; Gqi=x_qi ;
    x_qi=[1:qNum]./2.25 ;
    for qi=1:qNum
%         x_qi(qi)=1+(qi-1)*(99.999-1)/qNum ;
        Gqi(qi)=exp(-exp(-(x_qi(qi)-nu))/beta) ;
    end
    Gqi=[0.01;Gqi] ;
    clear nPoly q_Obs q_Sim
    nPoly=3 ;                                       % Polynomial order
    q_Obs(:,di)=quantile(Data_Obs,Gqi)' ;           % qNum quantiles of observation data
    q_Sim=quantile(Data_Sim,Gqi)' ;                 % qNum quantiles of simulated data
    Data_Sim_EGQM=zeros(length(Data_Sim),1) ;
    [Data_Sim_EGQM_qi,pos_Ref,p_EGQM]=DataClass_quantiles_Obs(qNum,nPoly,...
        q_Obs(:,di),Data_Obs,q_Sim,Data_Sim,0) ;
% % %     [Data_Sim_EGQM_qi,pos_Ref,p_EGQM]=DataClass_quantiles_Sim(qNum,nPoly,...
% % %         q_Obs(:,di),Data_Obs,q_Sim,Data_Sim,0) ;
    posTOT=0 ; DataTOT=0 ;
    for qi=1:length(q_Obs(:,di))+1
        Data_Sim_EGQM(pos_Ref{qi,1})=Data_Sim_EGQM_qi{qi,1} ;
            posTOT=posTOT+length(pos_Ref{qi,1}) ;
            DataTOT=DataTOT+length(Data_Sim_EGQM_qi{qi,1}) ;
    end
        %%% Double checking if all points of the dataset have been incorporated
        if posTOT~=length(Data_Sim_EGQM) || DataTOT~=length(Data_Sim_EGQM)
            disp('Simulation data classiffication process is not correct! Lengths of the vector before and after are not consistent \n')
        end

    %% Statistical evaluation and visualisation: Q-Q plot, CFD plot, bias via barplots, Taylor diagram
    % Bias computation
    Bias_Vec=([Data_Obs Data_Sim Data_Sim_Delta Data_Sim_EQM1 Data_Sim_EQM2 Data_Sim_EGQM]-Data_Obs) ;
    FigNum=100+di ;h3=figure(FigNum) ; set(h3,'units','normalized','outerposition',[0 0 1 1]) ;
    subplot(6,12,[2:5 14:17])
    for bci=2:length(Bias_Vec(1,:))
        [occ,bins,binspdf_full,pdf_full]=histfit_mod(Bias_Vec(:,bci),250,'kernel') ;
        bar(bins,occ/sum(occ)*100,1)
        hold on
        plot(binspdf_full,pdf_full/sum(occ)*100,'r-','LineWidth',2);
    end
    xlabel('Bias [m]')
    ylabel('Occurrence [%]')
    grid on
    legend('Raw simulation (ERA5)','','Delta','','EQM1','','EQM2','','EQGM','Location','Best')
    NormBias_Vec=sum(abs([Data_Obs Data_Sim Data_Sim_Delta Data_Sim_EQM1 Data_Sim_EQM2 Data_Sim_EGQM]-Data_Obs))./length(Data_Obs) ;
    % Conditional biases
% % %         %%% Case 1: 100 quantiles
% % %         nQuant1=[0.05:0.05:0.8 0.825 0.85 0.875 0.9:0.01:0.98 0.99:0.001:0.999] ;
% % %         clear Bias_q_Obs1 Data_Obs_cond1 Data_Sim_EQM1_cond1 Data_Sim_cond1 Data_Sim_EGQM_cond1 CondBias_Vec1
% % %         Bias_q_Obs1=quantile(Data_Obs,nQuant1)' ;
% % %         [Data_Obs_cond1,Data_Sim_cond1,Data_Sim_EQM2_cond1,Data_Sim_EGQM_cond1,CondBias_Vec1]=CondBias(nQuant1, ...
% % %             Bias_q_Obs1(1:end-1),Data_Obs,Data_Sim,Data_Sim_EQM2,Data_Sim_EGQM,0) ;
    %%% Case 2: 5 quantiles
    nQuant2=[0.25 0.5 0.75 0.90 0.99];%[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999] ;
    QuantileDef={'1^{st} quartile','2^{nd} quartile','3^{rd} quartile','75^{th}-90^{th} centiles',...
        '90^{th}-99^{th} centiles', '99^{th} centile'} ;
% % %     QuantileDef={'1^{st} decile','2^{nd} decile','3^{rd} decile','4^{th} decile','5^{th} decile',...
% % %             '6^{th} decile','7^{th} decile','8^{th} decile','9^{th} decile',...
% % %             '99^{th} centile', '99.9^{th} centile','rest'} ;
    clear Bias_q_Obs2 Data_Obs_cond2 Data_Sim_EQM2_cond2 Data_Sim_cond2 Data_Sim_EGQM_cond2 CondBias_Vec2
    Bias_q_Obs2=quantile(Data_Obs,nQuant2)' ;
    [Data_Obs_cond2,Data_Sim_cond2,Data_Sim_EQM2_cond2,Data_Sim_EGQM_cond2,CondBias_Vec2]=CondBias(nQuant2,QuantileDef, ...
        Bias_q_Obs2,Data_Obs,Data_Sim,Data_Sim_EQM2,Data_Sim_EGQM,FigNum) ;
    %%% 
    BarAlpha=[0.7 0.5 0.3] ;
    figure(FigNum) ;
    for ii=1:length(CondBias_Vec2)
        for jj=1:3
            clear occ bins binspdf_full pdf_full
            subplot(6,12,(3:5)*12+ii)
            hold on
            [occ,bins,binspdf_full,pdf_full]=histfit_mod(CondBias_Vec2{ii,1}(:,jj+1),250,'kernel') ;
%             h=bar(bins,occ/sum(occ)*100,1);
%                 h.FaceAlpha=BarAlpha(jj) ;
%                 hold on
%                 plot(binspdf,pdf/sum(occ)*100,'LineWidth',2);
            area(binspdf_full,pdf_full/sum(occ)*100,'FaceAlpha',BarAlpha(jj),'EdgeColor','k','LineStyle','-') ;
            title(QuantileDef{ii}) ;
            if ii==1 && jj==3
                legend('Raw simulations','EQM2','EGQM') ;
                xlabel('Bias [m]');
            end
            ylabel('PDF [-]');
            grid on
        end
        camroll(-270) ;
%         ylim([0 2.5])
%         xlim([0 5])
    end
    if SaveData
        saveas(h3,[ResultsFolderName,'/CondBias_',VarNames_var{1,di}(1:2)],'fig') ;
        saveas(h3,[ResultsFolderName,'/CondBias_',VarNames_var{1,di}(1:2)],'svg') ;
    end
    % Root mean Square Deviation (RMSD)
    RMSD_Vec=[0,...
        sqrt(mean((Data_Sim - Data_Obs).^2)),...
        sqrt(mean((Data_Sim_Delta - Data_Obs).^2)),...
        sqrt(mean((Data_Sim_EQM1 - Data_Obs).^2)),...
        sqrt(mean((Data_Sim_EQM2 - Data_Obs).^2)),...
        sqrt(mean((Data_Sim_EGQM - Data_Obs).^2))] ;
    % Standard Deviation (SD)
    SD_Vec=[std(Data_Obs),...
        std(Data_Sim),...
        std(Data_Sim_Delta),...
        std(Data_Sim_EQM1),...
        std(Data_Sim_EQM2),...
        std(Data_Sim_EGQM)] ;
    % Pearson Correlation (PC)
    PC_Vec=[corrcoef(Data_Obs,Data_Obs),...
        corrcoef(Data_Sim,Data_Obs),...
        corrcoef(Data_Sim_Delta,Data_Obs),...
        corrcoef(Data_Sim_EQM1,Data_Obs),...
        corrcoef(Data_Sim_EQM2,Data_Obs),...
        corrcoef(Data_Sim_EGQM,Data_Obs)] ;
    BCstat=table(NormBias_Vec',RMSD_Vec',SD_Vec',PC_Vec(1,2:2:end)','VariableNames',{'Bias','RMSD','SD','PC'}) ;

    % Q-Q plot
    h4=figure (200+di+1) ;
    QQplot_BC(Gqi(1:30)*100,Data_Obs,Data_Sim,Data_Sim_Delta,Data_Sim_EQM1,Data_Sim_EQM2,Data_Sim_EGQM,NormBias_Vec) ;
    FigTitle=['Q-Q plot for ',VarNames_var{1,di}(1:2)] ;
    sgtitle(FigTitle)
    if SaveData
        saveas(h4,[ResultsFolderName,'/QQplot_',VarNames_var{1,di}(1:2)],'fig') ;
        saveas(h4,[ResultsFolderName,'/QQplot_',VarNames_var{1,di}(1:2)],'svg') ;
    end
    % CFD plot
    h5=figure (300+di) ;
    CFDplot_BC(Data_Obs,Data_Sim,Data_Sim_Delta,Data_Sim_EQM1,Data_Sim_EQM2,Data_Sim_EGQM,NormBias_Vec) ;
    FigTitle=['CDF plot for ',VarNames_var{1,di}(1:2)] ;
    sgtitle(FigTitle)
    if SaveData
        saveas(h5,[ResultsFolderName,'/CDFplot_',VarNames_var{1,di}(1:2)],'fig') ;
        saveas(h5,[ResultsFolderName,'/CDFplot_',VarNames_var{1,di}(1:2)],'svg') ;
    end
    % Taylor diagram visualisation
    h6=figure (400+di) ;
    LabelCell={'Observation','ERA5_{raw}','ERA5_{BC-Delta}','ERA5_{BC-EQM1}','ERA5_{BC-EQM2}','ERA5_{BC-EGQM}'} ;
    taylordiag(SD_Vec',RMSD_Vec',PC_Vec(1,2:2:end)', ...
        'markerLabel',LabelCell, 'markerLegend', 'on', 'styleSTD', '-', 'colOBS','r', 'markerObs','*', ...
        'markerSize',10, 'tickRMS',0:0.5:2,'limSTD',2,'tickRMSangle',115, 'showlabelsRMS', 'on', ...
        'titleRMS','on', 'titleOBS','Observation') ;
    hold on
    FigTitle=['Taylor diagram for ',VarNames_var{1,di}(1:2)] ;
    sgtitle(FigTitle)
    if SaveData
        saveas(h6,[ResultsFolderName,'/TaylorDiagram_',VarNames_var{1,di}(1:2)],'fig') ;
        saveas(h6,[ResultsFolderName,'/TaylorDiagram_',VarNames_var{1,di}(1:2)],'svg') ;
    end
    %
    BCstats{di,1}=BCstat ;
    Data_Sim_corrected{di,1}=[Data_Obs,Data_Sim,Data_Sim_Delta,Data_Sim_EQM1,Data_Sim_EQM2,Data_Sim_EGQM] ;
    polyParams_EQM2{di,1}=p_EQM2 ;
    %
    a=table2array(MetoceanData_Sim(:,1)); b=table2array(MetoceanData_Obs(:,di+1)) ; a(isnan(b))=[] ;
    h7=figure(500+di); set(h7,'units','normalized','outerposition',[0 0 1 1]) ;
    plot(table2array(MetoceanData_Obs(:,1)),table2array(MetoceanData_Obs(:,di+1)),'-')
    hold on
    plot(table2array(MetoceanData_Sim(:,1)),table2array(MetoceanData_Sim(:,di+1)),'--')
    plot(a,Data_Sim_corrected{di,1}(:,end-1),'--')
    plot(a,Data_Sim_corrected{di,1}(:,end),'--')
    grid on
    datetick('x','yyyy')
    xlabel('Time [Years]')
    legend('Observation','Raw','BC-EQM2','BC-EGQM')
    if SaveData
        saveas(h7,[ResultsFolderName,'/CorrectedTimeSeries_',VarNames_var{1,di}(1:2)],'fig') ;
        saveas(h7,[ResultsFolderName,'/CorrectedTimeSeries_',VarNames_var{1,di}(1:2)],'svg') ;
    end
    %%%
    if di==1
        Hs=[datenum(Datetime_var) Data_Obs Data_Sim Data_Sim_EQM2 Data_Sim_EGQM] ;
    elseif di==2
        Tp=[datenum(Datetime_var) Data_Obs Data_Sim Data_Sim_EQM2 Data_Sim_EGQM] ;
    else
        Uw=[datenum(Datetime_var) Data_Obs Data_Sim Data_Sim_EQM2 Data_Sim_EGQM] ;
    end
% % %     Overall_1to99{di,1}=[[Data_Obs_cond2{1,1};Data_Obs_cond2{2,1};Data_Obs_cond2{3,1};Data_Obs_cond2{4,1};Data_Obs_cond2{5,1}],...
% % %         [Data_Sim_cond2{1,1};Data_Sim_cond2{2,1};Data_Sim_cond2{3,1};Data_Sim_cond2{4,1};Data_Sim_cond2{5,1}],...
% % %         [Data_Sim_EQM2_cond2{1,1};Data_Sim_EQM2_cond2{2,1};Data_Sim_EQM2_cond2{3,1};Data_Sim_EQM2_cond2{4,1};Data_Sim_EQM2_cond2{5,1}],...
% % %         [Data_Sim_EGQM_cond2{1,1};Data_Sim_EGQM_cond2{2,1};Data_Sim_EGQM_cond2{3,1};Data_Sim_EGQM_cond2{4,1};Data_Sim_EGQM_cond2{5,1}]] ;
    Overall{di,1}=[Data_Obs,Data_Sim,Data_Sim_EQM2,Data_Sim_EGQM] ;
    Extremes_99{di,1}=[[Data_Obs_cond2{end,1}],...
        [Data_Sim_cond2{end,1}],...
        [Data_Sim_EQM2_cond2{end,1}],...
        [Data_Sim_EGQM_cond2{end,1}]] ;
end

%% Assessment of the ORE resource
%%% PDF comarison: overall (0-99th percentiles) and extreme events (99th-100th perentiles) for Hs & Uw
PDFscore=zeros(4,4) ; DAV=zeros(4,2) ;
BarAlpha=[0.9 0.6 0.4 0.25] ; AreaColor= [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];
Overlapcolor={'r','y','m'} ; OverlapLineStyle={'--','-.',':'} ;
Nbins=250 ;
h8=figure(601) ;
ij=0 ;
for ii=1:2:3
    ij=ij+1 ;
    % Full dataset
    clear binspdf_full pdf_full PDFoverlap0_full PDFscore0_full DAV0_full
    [binspdf_full,pdf_full,PDFoverlap0_full,PDFscore0_full,DAV0_full]=SkillScore_fnc(Nbins,Overall{ii,1},max(max(Overall{ii,1}))+0.25*max(max(Overall{ii,1}))) ;
    subplot(2,2,ij)
    area(binspdf_full(:,1),pdf_full(:,1),'FaceAlpha',BarAlpha(1),'FaceColor',AreaColor(1,:),'EdgeColor','k','LineStyle','-')
    hold on 
    for jj=2:length(PDFscore0_full)
        area(binspdf_full(:,jj),pdf_full(:,jj),'FaceAlpha',BarAlpha(jj),'FaceColor',AreaColor(jj,:),'EdgeColor','k','LineStyle','-') ;
        plot(binspdf_full(:,jj),PDFoverlap0_full(:,jj),'Color',Overlapcolor{jj-1},'LineStyle',OverlapLineStyle{jj-1},'LineWidth',3)
        if ii==1 && jj==4
            title('Overall wave height distribution ')
            xlabel('H_s [m]')
            ylabel('PDF [-]')
            %%
        elseif ii==3 && jj==4
            title('Overall wind speed distribution ')
            xlabel('U_w [m/s]')
            ylabel('PDF [-]')
        end
    end
    hold off
    PDFscore(ij,:)=PDFscore0_full ; DAV(ij,:)=round(DAV0_full*100) ;
    legend('Observation',['Raw re-analysis (',num2str(PDFscore(ij,2)),' - N/A)'],...
        ['EQM (',num2str(PDFscore(ij,3)),' & ',num2str(DAV(ij,1)),'%)'],...
        ['EGQM (',num2str(PDFscore(ij,4)),' & ',num2str(DAV(ij,2)),'%)'])
    % 99th percentile
    clear binspdf_99 pdf_99 PDFoverlap0_99 PDFscore0_99 DAV0_99
    [binspdf_99,pdf_99,PDFoverlap0_99,PDFscore0_99,DAV0_99]=SkillScore_fnc(Nbins,Extremes_99{ii,1},max(max(Extremes_99{ii,1}))+0.25*max(max(Extremes_99{ii,1}))) ;
    subplot(2,2,ij+2)
    area(binspdf_99(:,1),pdf_99(:,1),'FaceAlpha',BarAlpha(1),'FaceColor',AreaColor(1,:),'EdgeColor','k','LineStyle','-')
    hold on 
    for jj=2:length(PDFscore0_99)
        area(binspdf_99(:,jj),pdf_99(:,jj),'FaceAlpha',BarAlpha(jj),'FaceColor',AreaColor(jj,:),'EdgeColor','k','LineStyle','-') ;
        plot(binspdf_99(:,jj),PDFoverlap0_99(:,jj),'Color',Overlapcolor{jj-1},'LineStyle',OverlapLineStyle{jj-1},'LineWidth',3)
        if ii==1 && jj==4
            title('Wave height extremes (99^{th} centile) distribution ')
            xlabel('H_s [m]')
            ylabel('PDF [-]')
            %%
        elseif ii==3 && jj==4
            title('Wind speed extremes (99^{th} centile) distribution ')
            xlabel('U_w [m/s]')
            ylabel('PDF [-]')
        end
    end
    hold off
    PDFscore(ij+2,:)=PDFscore0_99 ; DAV(ij+2,:)=round(DAV0_99*100) ;
    legend('Observation',['Raw re-analysis (',num2str(PDFscore(ij+2,2)),' - N/A)'],'',...
        ['EQM (',num2str(PDFscore(ij+2,3)),' & ',num2str(DAV(ij+2,1)),'%)'],'',...
        ['EGQM (',num2str(PDFscore(ij+2,4)),' & ',num2str(DAV(ij+2,2)),'%)'],'')
end
DAVmean=[(DAV(1,1)+DAV(3,1))/2 (DAV(1,2)+DAV(3,2))/2;(DAV(2,1)+DAV(4,1))/2 (DAV(2,2)+DAV(4,2))/2] ;
    set(h8,'units','normalized','outerposition',[0 0 1 1]) ;
    if SaveData
        saveas(h8,[ResultsFolderName,'/PDFscore'],'fig') ;
        saveas(h8,[ResultsFolderName,'/PDFscore'],'svg') ;
    end

%%% Power density (P)
if ~isempty(setdiff(Tp(:,1),Hs(:,1))) || ~isempty(setdiff(Hs(:,1),Tp(:,1)))
    [~,idx_Hs]=setdiff(Hs(:,1),Tp(:,1)) ;
        if ~isempty(idx_Hs)
            Hs(idx_Hs,:)= [] ;
        end
    [~,idx_Tp]=setdiff(Tp(:,1),Hs(:,1)) ;
        if ~isempty(idx_Tp)
            Tp(idx_Tp,:)= [] ;
        end
end
rho_H20=1025 ; rho_air=1.225 ;
alpha=0.9 ;
P_wave=zeros(size(Hs(:,1:end-1))) ; P_wind=zeros(size(Uw(:,1:end-1))) ;
for ii=2:length(Hs(1,:))
    % Wave power density (P_{wave})
    P_wave(:,ii-1)=0.49/alpha*Hs(:,ii).^2.*Tp(:,ii) ;
    % Wind power density (P_{wind})
    P_wind(:,ii-1)=0.5*rho_air*Uw(:,ii).^3 ;
end
P_wave_av=mean(P_wave) ;
P_wind_av=mean(P_wind) ;

%%% Coefficient of Variation (COV)
P_wave_sd=[std(P_wave(:,1)) std(P_wave(:,2)) std(P_wave(:,3)) std(P_wave(:,4))] ;
    COV_wave=P_wave_sd./P_wave_av ;
P_wind_sd=[std(P_wind(:,1)) std(P_wind(:,2)) std(P_wind(:,3)) std(P_wind(:,4))] ;
    COV_wind=P_wind_sd./P_wind_av ;
MarkerStyles={'o','s','d','h'} ;
h9=figure(602) ;
for ii=1:length(P_wave_av)
    subplot 121
    hold on
    h1=plot(COV_wave(ii),P_wave_av(ii),'Marker',MarkerStyles{ii},'MarkerSize',10);
    set(h1,'MarkerFaceColor',get(h1,'Color'));
    subplot 122
    hold on
    h2=plot(COV_wind(ii),P_wind_av(ii),'Marker',MarkerStyles{ii},'MarkerSize',10);
    set(h2,'MarkerFaceColor',get(h2,'Color'));
end
subplot 121
hold on
xlabel('CoV [-]')
ylabel('J_{wave} (kW/m)')
grid on
legend('Observation','Raw re-analysis','EQM','EGQM')
subplot 122
hold on
xlabel('CoV [-]')
ylabel('J_{wind} (kW/m^2)')
grid on
if SaveData
    saveas(h9,[ResultsFolderName,'/JvsCOV'],'fig') ;
    saveas(h9,[ResultsFolderName,'/JvsCOV'],'svg') ;
end

%%% Long-term trend: Monthly mean calculation & trend identification
% Hs
[mm_vec_Hs,Trend_Hs,slope_Hs,CI_lb_Hs,CI_ub_Hs,h10]=LongTermTrend(1,Hs,'H_s','m') ;
if SaveData
    saveas(h10,[ResultsFolderName,'/Trend_Hs'],'fig') ;
    saveas(h10,[ResultsFolderName,'/Trend_Hs'],'svg') ;
end
% Tp
[mm_vec_Tp,Trend_Tp,slope_Tp,CI_lb_Tp,CI_ub_Tp,h11]=LongTermTrend(2,Tp,'T_p','s') ;
if SaveData
    saveas(h11,[ResultsFolderName,'/Trend_Tp'],'fig') ;
    saveas(h11,[ResultsFolderName,'/Trend_Tp'],'svg') ;
end
% Uw
[mm_vec_Uw,Trend_Uw,slope_Uw,CI_lb_Uw,CI_ub_Uw,h12]=LongTermTrend(3,Uw,'U_w','m/s') ;
if SaveData
    saveas(h12,[ResultsFolderName,'/Trend_Uw'],'fig') ;
    saveas(h12,[ResultsFolderName,'/Trend_Uw'],'svg') ;
end
% Jwave
[mm_vec_Jwave,Trend_Jwave,slope_Jwave,CI_lb_Jwave,CI_ub_Jwave,h13]=LongTermTrend(4,[Hs(:,1) P_wave],'J_{wave}','kW/m') ;
if SaveData
    saveas(h13,[ResultsFolderName,'/Trend_Jwave'],'fig') ;
    saveas(h13,[ResultsFolderName,'/Trend_Jwave'],'svg') ;
end
% Jwind
[mm_vec_Jwind,Trend_Jwind,slope_Jwind,CI_lb_Jwind,CI_ub_Jwind,h14]=LongTermTrend(5,[Uw(:,1) P_wind],'J_{wind}','kW/m^2') ;
if SaveData
    saveas(h14,[ResultsFolderName,'/Trend_Jwind'],'fig') ;
    saveas(h14,[ResultsFolderName,'/Trend_Jwind'],'svg') ;
end

%% Save results
save([ResultsFolderName([1:8 17:end]),'_CorrectedData'],...
    'MetoceanData_Obs','MetoceanData_Sim','Hs','Tp','Uw','P_wave','P_wind',...
    'BCstats','PDFscore','DAV','DAVmean')