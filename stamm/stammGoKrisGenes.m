function terms=stammGoKrisGenes()
% STAMMGOKRISGENES List of interesting GO categories identified by Kris

% Cell signaling
% (11)JAK−STAT cascade
% (11)negative regulation of MAPKKK cascade
% (29)negative regulation of intracellular protein kinase cascade
% (36)MAPKKK cascade
% (33)regulation of MAPKKK cascade
% (55)regulation of small GTPase mediated signal transduction
% (50)regulation of cell division
terms=[7259; 43409; 10741; 165; 43408; 51056; 51302];

% cell migration/adhesion story
% (16)substrate−bound cell migration
% (5)heterophilic cell−cell adhesion
% (7)cell−matrix adhesion
% (316)cell surface receptor linked signaling pathway
% (0)integrin−mediated signaling pathway
terms=[terms; 6929; 7157; 7160; 7166; 7229];

% Proliferation and death story
% (357)cell cycle
% (100)cell death
% (76)apoptosis
% (3)anti−apoptosis
% (25)developmental cell growth
% (172)cell proliferation
% (32)induction of programmed cell death
terms=[terms; 7049; 8219; 6915; 6916; 48588; 8283; 12502];

% transdifferentation / nonpluripotent fates / niche formation?
% (294)embryonic morphogenesis
% (305)heart development
% (19)embryonic placenta development
% (213)regulation of leukocyte activation
% (20)axis elongation
% (0)patterning of blood vessels
% (59)epithelial tube formation
% (27)mammary gland epithelium development
% (34)muscle cell development
% (2)bone remodeling
% (460)embryonic development
% (444)generation of neurons
% (143)cell morphogenesis involved in differentiation
% (157)neuron development
% (27)skeletal muscle organ development
% (6)ectoderm development
% (235)immune system development
terms=[terms; 48598; 7507; 1892; 2694; 3401; 1569; 72175; 61180; 55001;
       46849; 9790; 48699; 902; 48666; 60538; 7398; 2520];

% transcriptional regulation
% (156)chromatin organization
% (25)negative regulation of transcription from RNA polymerase II promoter
% (85)transcription, DNA−dependent
% (129)chromatin modification
% (254)chromosome organization


terms=[terms; 6325; 122; 6351; 16568; 51276];
