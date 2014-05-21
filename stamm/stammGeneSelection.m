function genes=stammGeneSelection(x)
% STAMMGENESELECTION Return accession numbers of sets of core genes.

switch x
  case 1
    % Regulated genes from Fig 1 of Samavarchi-Tehrani 2010 + bonus genes
% $$$     genes={'Inadl'; 'Cldn11'; 'Cdh1'; 'Epcam'; 'Ocln'; 'Crb3'; 'Esrp1';
% $$$            'Cldn3'; 'Cldn4'; 'Cldn7'; 'Terf1'; 'Alpl'; 'Gdf3'; 'Tdgf1';
% $$$            'Sall4'; 'Utf1'; 'Zfp42'; 'Esrrb'; 'Lin28'; 'Prdm1'; 'Cdh2';
% $$$            'Ncam1'; 'Dnmt3l'; 'Zeb1'; 'Zeb2'; 'Thy1'; 'Dppa4';
% $$$            'Pecam1'; 'Snai2'; 'Tcf7'; 'Klf2'; 'Dazl'; 'Dppa2'; 'Dppa3';
% $$$            'Nodal'; 'Snai1'; 'Lefty1'; 'Gata4'; 'Tcl1'; 'Nanog' };
    genes={'NM_172696'; 'NM_008770'; 'NM_009864';
           'NM_008532'; 'NM_008756'; 'NM_177638'; 'NM_194055'; 
           'NM_009902'; 'NM_009903'; 'NM_016887'; 'NM_009352'; 
           'NM_007431'; 'NM_008108'; 'NM_011562'; 'NM_175303'; 
           'NM_009482'; 'NM_009556'; 'NM_001159500'; 'NM_145833'; 
           'NM_007548'; 'NM_007664'; 'NM_001113204'; 'NM_001081695'; 
           'NM_011546'; 'NM_015753'; 'NM_009382'; 'NM_028610'; 
           'NM_008816'; 'NM_011415'; 'NM_009331'; 'NM_008452'; 
           'NM_010021'; 'NM_028615'; 'NM_139218'; 'NM_013611'; 
           'NM_011427'; 'NM_010094'; 'NM_008092'; 'NM_009337'; 
           'NM_028016'};
  case 2
    % as above plus, known iPSC blockers and key genes
% $$$     genes={'Inadl'; 'Cldn11'; 'Cdh1'; 'Epcam'; 'Ocln'; 'Crb3'; 'Esrp1';
% $$$            'Cldn3'; 'Cldn4'; 'Cldn7'; 'Terf1'; 'Alpl'; 'Gdf3'; 'Tdgf1';
% $$$            'Sall4'; 'Utf1'; 'Zfp42'; 'Esrrb'; 'Lin28'; 'Prdm1'; 'Cdh2';
% $$$            'Ncam1'; 'Dnmt3l'; 'Zeb1'; 'Zeb2'; 'Thy1'; 'Dppa4';
% $$$            'Pecam1'; 'Snai2'; 'Tcf7'; 'Klf2'; 'Dazl'; 'Dppa2'; 'Dppa3';
% $$$            'Nodal'; 'Snai1'; 'Lefty1'; 'Gata4'; 'Tcl1'; 'Nanog';
% $$$           'Tgfb1'; 'Tgfb2'; 'Cdkn2a'; 'Cdkn2b'; 'Dppa5a'};
    genes={'NM_172696'; 'NM_008770'; 'NM_009864';
           'NM_008532'; 'NM_008756'; 'NM_177638'; 'NM_194055'; 
           'NM_009902'; 'NM_009903'; 'NM_016887'; 'NM_009352'; 
           'NM_007431'; 'NM_008108'; 'NM_011562'; 'NM_175303'; 
           'NM_009482'; 'NM_009556'; 'NM_001159500'; 'NM_145833'; 
           'NM_007548'; 'NM_007664'; 'NM_001113204'; 'NM_001081695'; 
           'NM_011546'; 'NM_015753'; 'NM_009382'; 'NM_028610'; 
           'NM_008816'; 'NM_011415'; 'NM_009331'; 'NM_008452'; 
           'NM_010021'; 'NM_028615'; 'NM_139218'; 'NM_013611'; 
           'NM_011427'; 'NM_010094'; 'NM_008092'; 'NM_009337'; 
           'NM_028016'; 'NM_011577'; 'NM_009367'; 'NM_009877'; 
           'NM_007670'; 'NM_025274'};
  case 3
    % as above plus Oct4, Sox2, Klf4, c-Myc
% $$$     genes={'Inadl'; 'Cldn11'; 'Cdh1'; 'Epcam'; 'Ocln'; 'Crb3'; 'Esrp1';
% $$$            'Cldn3'; 'Cldn4'; 'Cldn7'; 'Terf1'; 'Alpl'; 'Gdf3'; 'Tdgf1';
% $$$            'Sall4'; 'Utf1'; 'Zfp42'; 'Esrrb'; 'Lin28'; 'Prdm1'; 'Cdh2';
% $$$            'Ncam1'; 'Dnmt3l'; 'Zeb1'; 'Zeb2'; 'Thy1'; 'Dppa4';
% $$$            'Pecam1'; 'Snai2'; 'Tcf7'; 'Klf2'; 'Dazl'; 'Dppa2'; 'Dppa3';
% $$$            'Nodal'; 'Snai1'; 'Lefty1'; 'Gata4'; 'Tcl1'; 'Nanog';
% $$$            'Tgfb1'; 'Tgfb2'; 'Cdkn2a'; 'Cdkn2b'; 'Dppa5a'; 'Pou5f1';
% $$$            'Sox2'; 'Klf4'; 'Myc'};
    genes={'NM_172696'; 'NM_008770'; 'NM_009864';
           'NM_008532'; 'NM_008756'; 'NM_177638'; 'NM_194055'; 
           'NM_009902'; 'NM_009903'; 'NM_016887'; 'NM_009352'; 
           'NM_007431'; 'NM_008108'; 'NM_011562'; 'NM_175303'; 
           'NM_009482'; 'NM_009556'; 'NM_001159500'; 'NM_145833'; 
           'NM_007548'; 'NM_007664'; 'NM_001113204'; 'NM_001081695'; 
           'NM_011546'; 'NM_015753'; 'NM_009382'; 'NM_028610'; 
           'NM_008816'; 'NM_011415'; 'NM_009331'; 'NM_008452'; 
           'NM_010021'; 'NM_028615'; 'NM_139218'; 'NM_013611'; 
           'NM_011427'; 'NM_010094'; 'NM_008092'; 'NM_009337'; 
           'NM_028016'; 'NM_011577'; 'NM_009367'; 'NM_009877'; 
           'NM_007670'; 'NM_025274'; 'NM_013633'; 'NM_011443'; 
           'NM_010637'; 'NM_010849'};
    
  case 4
    % as above plus matching genes from Fluidigm data
% $$$     genes={'Inadl'; 'Cldn11'; 'Cdh1'; 'Epcam'; 'Ocln'; 'Crb3'; 'Esrp1';
% $$$            'Cldn3'; 'Cldn4'; 'Cldn7'; 'Terf1'; 'Alpl'; 'Gdf3'; 'Tdgf1';
% $$$            'Sall4'; 'Utf1'; 'Zfp42'; 'Esrrb'; 'Lin28'; 'Prdm1'; 'Cdh2';
% $$$            'Ncam1'; 'Dnmt3l'; 'Zeb1'; 'Zeb2'; 'Thy1'; 'Dppa4';
% $$$            'Pecam1'; 'Snai2'; 'Tcf7'; 'Klf2'; 'Dazl'; 'Dppa2'; 'Dppa3';
% $$$            'Nodal'; 'Snai1'; 'Lefty1'; 'Gata4'; 'Tcl1'; 'Nanog';
% $$$            'Tgfb1'; 'Tgfb2'; 'Cdkn2a'; 'Cdkn2b'; 'Dppa5a'; 'Pou5f1';
% $$$            'Sox2'; 'Klf4'; 'Myc'; 'Col5a2'; 'Fut4'; 'Tbx3'; 'Stat3';
% $$$            'Ctnnbl1'; 'Lifr'; 'Slc2a1'; 'Nes'; 'Bmi1'; 'Ctcf';
% $$$            'Myst4'; 'Nr6a1'; 'Dnmt3b'; 'Ezh2'; 'Gsk3b'; 'Csnk2a1';
% $$$            'Kdm1'; 'Dnmt1'; 'Prmt7'; 'Cdc20'; 'Mad2l1'; 'Ccnf'; 'Fgf5';
% $$$            'Fgf4'; 'Jag1'; 'Notch1'; 'Bub1'; 'Hprt1'};
    genes={'NM_172696'; 'NM_008770'; 'NM_009864';
           'NM_008532'; 'NM_008756'; 'NM_177638'; 'NM_194055'; 
           'NM_009902'; 'NM_009903'; 'NM_016887'; 'NM_009352'; 
           'NM_007431'; 'NM_008108'; 'NM_011562'; 'NM_175303'; 
           'NM_009482'; 'NM_009556'; 'NM_001159500'; 'NM_145833'; 
           'NM_007548'; 'NM_007664'; 'NM_001113204'; 'NM_001081695'; 
           'NM_011546'; 'NM_015753'; 'NM_009382'; 'NM_028610'; 
           'NM_008816'; 'NM_011415'; 'NM_009331'; 'NM_008452'; 
           'NM_010021'; 'NM_028615'; 'NM_139218'; 'NM_013611'; 
           'NM_011427'; 'NM_010094'; 'NM_008092'; 'NM_009337'; 
           'NM_028016'; 'NM_011577'; 'NM_009367'; 'NM_009877'; 
           'NM_007670'; 'NM_025274'; 'NM_013633'; 'NM_011443'; 
           'NM_010637'; 'NM_010849'; 'NM_007737'; 'NM_010242';
           'NM_011535'; 'NM_213659'; 'NM_025680'; 'NM_013584';
           'NM_011400'; 'NM_016701'; 'NM_007552'; 'NM_181322';
           'NM_017479'; 'NM_010264'; 'NM_001003961'; 'NM_007971';
           'NM_019827'; 'NM_007788'; 'NM_133872'; 'NM_010066';
           'NM_145404'; 'NM_023223'; 'NM_019499'; 'NM_007634';
           'NM_010203'; 'NM_010202'; 'NM_013822'; 'NM_008714';
           'NM_001113179'; 'NM_013556'};
  case 5
    % as above with Oct4, Sox2, c-Myx, Klf4 removed and with
    % Ccnd1, Trim28, Nr0b1, Zfp248, Suz12, Eed, Phc1, Rnf2, Tcf3,
    % Klf5, Rest
    genes={'NM_172696'; 'NM_008770'; 'NM_009864';
           'NM_008532'; 'NM_008756'; 'NM_177638'; 'NM_194055'; 
           'NM_009902'; 'NM_009903'; 'NM_016887'; 'NM_009352'; 
           'NM_007431'; 'NM_008108'; 'NM_011562'; 'NM_175303'; 
           'NM_009482'; 'NM_009556'; 'NM_001159500'; 'NM_145833'; 
           'NM_007548'; 'NM_007664'; 'NM_001113204'; 'NM_001081695'; 
           'NM_011546'; 'NM_015753'; 'NM_009382'; 'NM_028610'; 
           'NM_008816'; 'NM_011415'; 'NM_009331'; 'NM_008452'; 
           'NM_010021'; 'NM_028615'; 'NM_139218'; 'NM_013611'; 
           'NM_011427'; 'NM_010094'; 'NM_008092'; 'NM_009337'; 
           'NM_028016'; 'NM_011577'; 'NM_009367'; 'NM_009877'; 
           'NM_007670'; 'NM_025274'; 'NM_007737'; 'NM_010242';
           'NM_011535'; 'NM_213659'; 'NM_025680'; 'NM_013584';
           'NM_011400'; 'NM_016701'; 'NM_007552'; 'NM_181322';
           'NM_017479'; 'NM_010264'; 'NM_001003961'; 'NM_007971';
           'NM_019827'; 'NM_007788'; 'NM_133872'; 'NM_010066';
           'NM_145404'; 'NM_023223'; 'NM_019499'; 'NM_007634';
           'NM_010203'; 'NM_010202'; 'NM_013822'; 'NM_008714';
           'NM_001113179'; 'NM_013556'; 'NM_007631'; 'NM_011588';
           'NM_007430'; 'NM_025788'; 'NM_028335'; 'NM_199196';
           'NM_021876'; 'NM_007905'; 'NM_011277'; 'NM_001079822';
           'NM_009769'; 'NM_011263'
          };
  otherwise
    error('Unknown gene selection set');
end
