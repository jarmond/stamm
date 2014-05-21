function terms=stammGoMarioGenes(shorter)
% STAMMGOMARIOGENES List of interesting categories identified by Mario.

if nargin<1
    shorter=0;
end

% cell migration/adhesion story
% $$$ locomotion
% $$$ cell motility
% $$$ localization of cell
% $$$ epithelium development
% $$$ morphogenesis of an epithelium
terms=[40011; 48870; 51674; 60429; 2009];

% transdifferentation / nonpluripotent fates / niche formation?
% $$$ embryonic morphogenesis
% $$$ tissue morphogenesis
% $$$ embryonic development
% $$$ immune system development
% $$$ heart development
% $$$ post−embryonic development
% $$$ cell fate commitment
% $$$ generation of neurons
% $$$ organ morphogenesis
% $$$ cell development
% $$$ negative regulation of developmental process
% $$$ regulation of anatomical structure morphogenesis
% $$$ sensory organ development
% $$$ neurogenesis
% $$$ anatomical structure formation involved in morphogenesis
terms=[terms; 48598; 48729; 9790; 2520; 7507; 9791; 45165;
       48699; 9887; 48468; 51093; 22603; 7423; 22008; 48646];

% Cell signalling
% $$$ cell surface receptor linked signaling pathway
% $$$ regulation of cell communication
% $$$ negative regulation of signal transduction
% $$$ negative regulation of signaling process
% $$$ response to endogenous stimulus
% $$$ response to abiotic stimulus
% $$$ negative regulation of response to stimulus
% $$$ negative regulation of cell communication
terms=[terms; 7166; 10646; 9968; 23057; 9719; 9628; 48585; 10648];

% Proliferation and death story
% $$$ cell cycle process
% $$$ regulation of cell cycle process
% $$$ cell cycle
% $$$ regulation of cell cycle
terms=[terms; 22402; 10564; 7049; 51726];
       
% transcriptional regulation
% $$$ regulation of transcription, DNA−dependent
% $$$ chromosome organization
terms=[terms; 6355; 51276];

% other
% $$$ macromolecule catabolic process
% $$$ organic acid biosynthetic process
% $$$ small molecule biosynthetic process
% $$$ carboxylic acid biosynthetic process
% $$$ regulation of phosphate metabolic process
% $$$ regulation of catalytic activity
% $$$ regulation of phosphorus metabolic process
% $$$ positive regulation of multicellular organismal process
% $$$ regulation of cytokine production
% $$$ regulation of leukocyte activation
% $$$ protein localization
% $$$ negative regulation of biosynthetic process
% $$$ negative regulation of metabolic process
% $$$ regulation of catabolic process
% $$$ negative regulation of macromolecule metabolic process
% $$$ cellular membrane organization
% $$$ cellular component assembly
% $$$ negative regulation of cellular metabolic process
% $$$ negative regulation of cellular biosynthetic process
% $$$ regulation of cellular protein metabolic process
% $$$ regulation of system process
% $$$ cellular component biogenesis
% $$$ intracellular transport
% $$$ regulation of immune response
% $$$ positive regulation of transport
% $$$ regulation of cellular component organization
% $$$ regulation of RNA metabolic process
% $$$ membrane organization
% $$$ regionalization
% $$$ monocarboxylic acid metabolic process
% $$$ establishment of protein localization
% $$$ organelle organization
% $$$ regulation of transport
% $$$ regulation of molecular function
% $$$ response to organic substance
% $$$ macromolecular complex assembly
% $$$ monosaccharide metabolic process
% $$$ cellular component morphogenesis
% $$$ macromolecular complex subunit organization
% $$$ macromolecule localization
% $$$ nucleobase, nucleoside and nucleotide metabolic process
% $$$ positive regulation of molecular function
% $$$ cellular protein localization
% $$$ negative regulation of multicellular organismal process
% $$$ cellular macromolecule localization
% $$$ regulation of cellular catabolic process
% $$$ reproductive cellular process
% $$$ cellular localization
% $$$ regulation of organelle organization
% $$$ small molecule catabolic process
% $$$ establishment of localization in cell
% $$$ regulation of cellular localization
% $$$ cell activation
% $$$ nucleoside phosphate metabolic process
% $$$ nucleotide metabolic process
% $$$ negative regulation of transport
if shorter==0
    terms=[terms; 9057; 16053; 44283; 46394; 19220; 50790; 51174;
       51240; 1817; 2694; 8104; 9890; 9892; 9894; 10605; 16044;
       22607; 31324; 31327; 32268; 44057; 44085; 46907; 50776;
       51050; 51128; 51252; 61024; 3002; 32787; 45184; 6996;
       51049; 65009; 10033; 65003; 5996; 32989; 43933; 33036;
       55086; 44093; 34613; 51241; 70727; 31329; 51641; 33043;
       44282; 51234; 60341; 1775; 6753; 9117; 51051];
end
