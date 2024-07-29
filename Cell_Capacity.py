import numpy as np
import pandas as pd
import time as timer
import seaborn as sns
import matplotlib.pyplot as plt

colors = ['#1F77B4', '#FF7F0E', '#2CA02C', "#D62728"]

Cm = ["chloramphenicol","Chloramphenicol"]
Tet = ["tetracycline efflux","Tetracycline efflux","TetA","Tet(A)","tetA","tetracycline-inactivating"]
MLS = ["macrolide","lincosamide","streptogramin"]
multidrug = ["Multidrug resistance","multidrug resistance","antibiotic resistance"]
blactam = ["lactamase","LACTAMASE","beta-lactam","oxacillinase","carbenicillinase","betalactam\\S*"]
glycopeptide = ["glycopeptide resistance","VanZ","vancomycin resistance","VanA","VanY","VanX",
              "VanH","streptothricin N-acetyltransferase"]
polypeptide = ["bacitracin","polymyxin B","phosphoethanolamine transferase",
             "phosphoethanolamine--lipid A transferase"]
Trim = ["trimethoprim","dihydrofolate reductase","dihydropteroate synthase"]
Sulf = ["sulfonamide","Sul1","sul1","sulphonamide"]
Qnr = ["quinolone","Quinolone","oxacin","qnr","Qnr"]
Kan = ["Aminoglycoside","aminoglycoside","streptomycin","Streptomycin","kanamycin","Kanamycin","tobramycin","Tobramycin","gentamicin","Gentamicin","neomycin","Neomycin",
     "16S rRNA (guanine(1405)-N(7))-methyltransferase","23S rRNA (adenine(2058)-N(6))-methyltransferase","spectinomycin 9-O-adenylyltransferase",
     "Spectinomycin 9-O-adenylyltransferase","Rmt"]
macrolide = ["macrolide","ketolide","Azithromycin","azithromycin","Clarithromycin","clarithromycin","Erythromycin","erythromycin","Erm","EmtA"]
antimicrobial = ["QacE","Quaternary ammonium","quaternary ammonium",
               "Quarternary ammonium","quartenary ammonium","fosfomycin","ribosomal protection",
               "rifampin ADP-ribosyl","azole resistance","antimicrob\\S*"]

Resistance = [blactam, Cm,Tet, macrolide, Qnr, Sulf, Trim, Kan, multidrug]
ResistType = [r"$\beta$-lactam", "Cm", "Tet", "Macrolide", "Qnr", "Sulf", "Trim", "AMG", "Multi-drug Resistant"]

try:
    df=pd.read_csv("Hawkey2022_cell_capacity.csv")
except:
    df_length=pd.read_csv("Hawkey2022_replicon_lengths.csv")
    df_cn=pd.read_csv("Hawkey2022_chromosome_plasmid_copy_numbers.csv")
    df_ar=pd.read_csv("Hawkey2022_ARG_copy_numbers.csv")
    
    cell_list=[]
    annotations=df_cn["AnnotationAccession"].values.tolist()
    for aa in annotations:
        if aa not in cell_list:
            cell_list.append(aa)
    df={"AnnotationAccession":[],"total BP":[],"MDR":[],"Resistome":[],"total CN":[],
        "weighted size":[],"average size":[],"normalized CN":[]}
    for j in range(len(ResistType)):
        df[ResistType[j]]=[]
    t0=timer.perf_counter()
    print("Search Started")
    for cell in cell_list:
        df1=df_cn[df_cn["AnnotationAccession"]==cell]
        plasmid_id_cn=df1.index[df1["SeqType"]=="plasmid"]
        plasmid_id_cn=np.array(plasmid_id_cn)

        BP=0
        CN=0
        Size=[]
        hot_vector = np.zeros(len(ResistType))
        for i in plasmid_id_cn:
            seq_ID=df_cn.loc[i,"SeqID"]
            cn=df_cn.loc[i,"CopyNumber"]
            if cn<1:
                continue
            len_ids=np.where(df_length["SeqID"]==seq_ID)[0]
            if len(len_ids)>1:
                print("found more than 1 replicon length, seq ID: ",seq_ID)
            length=df_length.loc[len_ids[0],"replicon_length"]
            BP+=cn*length
            CN+=cn
            Size.append(length)
            ar_ids=np.where(df_ar["SeqID"]==seq_ID)[0]

            for id in ar_ids:
                description=df_ar.loc[id,"product"]
                found = False
                for j in range(len(ResistType)):
                    for k in Resistance[j]:
                        if k in description:
                            hot_vector[j]=1
                            break
        if BP<=1:
            continue
        if CN==0:
            continue
        df["AnnotationAccession"] += [cell, ]
        df["total BP"]+=[BP,]
        df["total CN"]+=[CN,]
        df["weighted size"]+=[BP/CN,]
        avg_sz=10**np.mean(np.log10(Size))
        df["average size"]+=[avg_sz,]
        df["normalized CN"]+=[BP/avg_sz,]
        for j in range(len(ResistType)):
            df[ResistType[j]]+=[True if hot_vector[j]==1 else False]
        if np.sum(hot_vector)>2:
            df["MDR"]+=[True,]
        else:
            df["MDR"] += [False, ]
        if np.sum(hot_vector)>2:
            df["Resistome"]+=["MDR",]
        elif np.sum(hot_vector)==1:
            i=int(np.where(hot_vector==1)[0][0])
            df["Resistome"] += [ResistType[i], ]
        else:
            df["Resistome"] += ["Sus", ]
    df=pd.DataFrame(df)
    df.to_csv("Hawkey2022_cell_capacity.csv")


    
ho=["Sus","MDR","Trim",r"$\beta$-lactam","Tet","AMG"]
fig,axes=plt.subplots(1,2,figsize=(6,3))
ax1,ax2=axes
df0=pd.read_csv("Hawkey2022_plasmid_capacity.csv")
bw=5e4
sns.histplot(data=df0,x="plasmid capacity",hue="Resistome",multiple="stack",ax=ax1,
             linewidth=0.3,log_scale=True,hue_order=ho)
sns.histplot(data=df,x="total BP",hue="Resistome",multiple="stack",ax=ax2,
             linewidth=0.3,log_scale=True,hue_order=ho)#binwidth=1e5,
ax1.set_xlim([1e3,5e6])
ax2.set_xlim([1e3,5e6])
ax1.legend().remove()
ax1.set_xlabel("plasmid bp")
ax2.set_xlabel("cell x-chromosomal bp")
fig.tight_layout()
fig.savefig("cell_capacity.png",dpi=220)

titles=["constant copy number","average size"]
xlabels=["weighted size","average size"]
ylabels=["total CN","normalized CN"]
zlabel = "total BP"
R388_size=np.array([33926,33926])
R388_cn=np.array([3.8,1.8])
pNUK73_size=np.array([5130,5130])
pNUK73_cn=np.array([11,9.56])
p_size=[R388_size,pNUK73_size]
p_cn=[R388_cn,pNUK73_cn]
for i in range(2):
    fig,axes=plt.subplots(1,2,figsize=(6,3))
    title=titles[i]
    xlabel=xlabels[i]
    ylabel=ylabels[i]
    fig.suptitle(title)
    ax1,ax2=axes
    sns.scatterplot(data=df,x=xlabel,y=ylabel,hue="Resistome",ax=ax1)
    xlim=ax1.get_xlim()
    for i in range(2):
        size = p_size[i]
        cn = p_cn[i]
        ax1.scatter(size, cn, facecolor=colors[2], edgecolor="k", linewidth=1)
        ax1.annotate("", xy=[size[1], cn[1]], xytext=[size[0], cn[0]],
                     arrowprops=dict(arrowstyle="->", color=colors[2]))
    x=np.array(df[xlabel])
    y=np.array(df[ylabel])
    z=np.polyfit(np.log10(x),np.log10(y),1)
    logxx=np.linspace(2.8,6)
    f=np.poly1d(z)
    logyy=f(logxx)
    xx=np.power(10,logxx)
    yy=np.power(10,logyy)
    ax1.plot(xx,yy,linewidth=2,c="k",zorder=10)
    ax1.text(x=5e4,y=3e1,s="y$\propto$x$^{%1.2f}$"%(z[0]))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_ylim([1,2e3])
    ax1.set_xlim([1e3,xlim[1]])
    ax1.legend(bbox_to_anchor=(0,1),loc="lower left",ncol=6,frameon=False,handlelength=0.2)#,title="multi-resistance")

    x=np.array(df[xlabel])
    y=np.array(df[zlabel])
    z=np.polyfit(np.log10(x),np.log10(y),1)
    logxx=np.linspace(2.8,6)
    f=np.poly1d(z)
    logyy=f(logxx)
    xx=np.power(10,logxx)
    yy=np.power(10,logyy)
    ax2.plot(xx,yy,linewidth=2,c="k",zorder=10)
    ax2.text(x=5e4,y=4e4,s="y$\propto$x$^{%1.2f}$"%(z[0]))
    sns.scatterplot(data=df,x=xlabel,y=zlabel,hue="Resistome",ax=ax2,legend=False)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim([1e3,xlim[1]])
    ax2.set_ylim([5e3,5e6])
    fig.subplots_adjust(wspace=0.3,bottom=0.22,top=0.8)
    fig.savefig("%s.png"%title,dpi=220)
plt.show()
