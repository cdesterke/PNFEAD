list.files()

data<-read.table("PNF_EAD.csv",sep="\t",h=T,na.string="NA")


library(Hmisc)

res<-describe(data)


y<-colnames(data)

equation<-paste(y,collapse="+")

data$outcome_PNF_EAD<-as.factor(data$outcome_PNF_EAD)

library(Publish)

data<-read.table("PNF_EAD.csv",sep="\t",h=T,na.string="NA")


u<-univariateTable(outcome_PNF_EAD~AGE+SEXE+INFARCTUS_MYOCARDE+INSUFFISANCE_CARDIAQUE+
PATHOLOGIE_VASCULAIRE+AVC_AIT+PATHO_PULMONAIRE_CHRONIQUE+DIABETE+
COMPLICATION_DIABETE+CANCER+VIH+ABO+Rhesus+INDICATION_TH+
SI_CIRRHOSE_INDICATION+CAUSE_CIRRHOSE+ENDROIT+GRADE_ACLF+
CREAT_microM+DIALYSE+IOT+INOTROPES+MELD_inscription+INR+BILI+
PLQ+ASCITE+ILA+HYDROTHORAX+THROMBOSE_PORTE+SHP+HTAP+EH+GRADE_EH+
TIPSS+ATCD_CHC+CHC_ACTIF_A_TH+AFP+AGE_DONNEUR+SEXE_DONNEUR+
POIDS_DONNEUR+TAILLE_DONNEUR+IMC_DONNEUR+CAUSE_MORT_ENCEPHALIQUE.DONNEUR+
ARRET_CARDIORESP_DONNEUR+NOMBRE_JOUR_REA_DONNEUR+CREAT_DONNEUR+
ASAT_DONNEUR+ALAT_DONNEUR+GGT_DONNEUR+PAL_DONNEUR+BILI_DONNEUR+
LACTATES_DONNEUR+NA_MGH_DONNEUR+PROXIMITE_PRELEVEMENT+TYPE_ORGANE+
INTEGRITE_FOIE+GREFFE_COMBINE+LIGAMENT_ARQUE+ANASTOMOSE_PORTE+
FERMETURE_SHUNT+THOMBECTOMIE+ANASTOMOSE_SUSHEPATIQUE+ANASTOMOSE_BILIAIRE+
TEMPS_ISCHEMIE_FROIDE+TEMPS_ISCHEMIE_TIEDE+TEMPS_ISCHEMIE_CHAUDE+
TEMPS_INTERVENTION+DIFFICULTE+CGR+CP+PFC+PERTES_SANGUINES+
SYNDROME_REPERFUSION+POIDS_GREFFON+POIDS_EXPLANT+CEC+
RECONTRUCTION_ARTERE+RECONTRUCTION_PORTE+CLAMPAGE_TOTAL_CAVE+VAC+
NBRE_ESP_PORTE_J0+J0_STEATOSE_MICRO+J0_STEATOSE_MACRO+J0_FIBROSE_PORTALE+
J0_FIBROSE_PERISINUSOIDALE+J0_FIBROSE_SIN_SYSTEMATISEE+J0_NECROSE+
J0_INFLAMMATION+J0_FER+J0_ISCHEMIE_REPERFUSION+NBRE_ESP_PORTE_REPERFUSION+
REPERFUSION_STEATOSE_MICRO+REPERFUSION_STEATOSE_MACRO+REPERFUSION_FIBROSE_PORTALE+
REPERFUSION_FIBROSE_PERISINUSOIDALE+REPERFUSION_FIBROSE_SIN_SYSTEMATISEE+
REPERFUSION_NECROSE_LOBULAIRE+REPERFUSION_NECROSE_VEINEUSE_PONCTUEE+
REPERFUSION_INFLAMMATION+REPERFUSION_FER+REPERFUSION_ISCHEMIE_REPERFUSION+
THROMBOSE_ARTERE+STENOSE_ARTERE+EVENEMENT_ARTERIEL+BILIsup100_J7+
CHOLANGITE_ISCHEMIQUE+THROMBOSE._SUSHEP+THROMBOSE_PORTE_POSTTH+STENOSE_PORTE+
STENOSE_BILIAIRE+FUITE_BILIAIRE+EVENEMENT_BILIAIRE+REPRISE_BLOC+REJET+DECES+
DECES_LIE_A_LA_GREFFE+GREFFON_FONCTIONNEL+NOUVELLE_TH+LACTATES_J0+NA_MGH_J1+
Na_J1+LACTATES_J1+FIO2_J1+PO2_J1+TP_J1+INR_J1+BILI_J1+CREAT_J1+ASAT_J1+ALAT_J1+
GGT_J1+PAL_J1+Hb_J1+GB_J1+PLQ_J1+DIALYSE_J1+IOT_J1+NA_MGH_J4+Na_J4+LACTATES_J4+
FIO2_J4+PO2_J4+TP_J4+INR_J4+BILI_J4+CREAT_J4+ASAT_J4+ALAT_J4+GGT_J4+PAL_J4+Hb_J4+
GB_J4+PLQ_J4+DIALYSE_J4+IOT_J4+NA_MGH_J7+Na_J7+LACTATES_J7+FIO2_J7+PO2_J7+TP_J7+
INR_J7+BILI_J7+CREAT_J7+ASAT_J7+ALAT_J7+GGT_J7+PAL_J7+Hb_J7+GB_J7+PLQ_J7+
ECHO_J7_AR+DIALYSE_J7+IOT_J7+NA_MGH_J14+Na_J14+LACTATES_J14+FIO2_J14+PO2_J14+
TP_J14+INR_J14+BILI_J14+CREAT_J14+ASAT_J14+ALAT_J14+GGT_J14+PAL_J14+Hb_J14+
GB_J14+PLQ_J14+DIALYSE_J14+IOT_J14+NA_MGH_J21+Na_J21+LACTATES_J21+FIO2_J21+
PO2_J21+TP_J21+INR_J21+BILI_J21+CREAT_J21+ASAT_J21+ALAT_J21+GGT_J21+PAL_J21+
Hb_J21+GB_J21+PLQ_J21+DIALYSE_J21+IOT_J21,data=data)


res<-summary(u)

res

write.table(res,file="PNF_EAD_univariable_results.tsv",sep="\t",row.names=F)

