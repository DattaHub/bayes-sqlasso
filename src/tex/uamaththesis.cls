\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{uamaththesis}[2014/9/28 Class for University of
    Arkansas math theses]
\newif\ifchapters
\chaptersfalse
\DeclareOption{chapters}{\chapterstrue}
\DeclareOption*{
  \ifchapters % Set by the class option.
    \PassOptionsToClass{\CurrentOption}{report}%
  \else
    \PassOptionsToClass{\CurrentOption}{article}%
  \fi
}
\ProcessOptions
\ifchapters
  \LoadClass[12pt]{report}
  \def\@makechapterhead#1{%
    {\parindent \z@ \raggedright \normalfont
      \ifnum \c@secnumdepth >\m@ne
          \normalfont\bfseries \@chapapp\space \thechapter
          \par\nobreak
          %\vskip 20\p@
      \fi
      \interlinepenalty\@M
      \normalfont \bfseries #1\par\nobreak
      %\vskip 40\p@
    }%
  }
  \def\@makeschapterhead#1{%
    {\parindent \z@ \raggedright
      \normalfont
      \interlinepenalty\@M
      \normalfont \bfseries  #1\par\nobreak
      %\vskip 40\p@
    }%
  }
\else
  \LoadClass[12pt]{article}
\fi

\renewcommand{\contentsname}{Table of Contents}
\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 0.1em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}
\providecommand\refname{}
% They seem to allow "Bibliography" now, not that it matters much.
\renewcommand{\refname}{References}
%
% Footnotes the same size as regular text
\long\def\@footnotetext#1{\insert\footins{%
    \reset@font\normalsize
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
    \protected@edef\@currentlabel{%
       \csname p@footnote\endcsname\@thefnmark
    }%
    \color@begingroup
      \@makefntext{\normalsize\selectfont%
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox}%
    \color@endgroup}}%

% Section headings same size as normal text
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   %{2.3ex \@plus.2ex}%
                                   {0.1ex \@plus.2ex}%
                                   {\normalfont\normalsize\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     %{1.5ex \@plus .2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\normalfont\normalsize\scshape}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\itshape}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\scshape}}

\newif\if@mainmatter
\newcommand\frontmatter{%
  \clearpage
  \@mainmatterfalse
  \pagestyle{empty}}
\newcommand\mainmatter{%
  \clearpage
  \@mainmattertrue
  \pagenumbering{arabic}%
  \pagestyle{plain}}
\newcommand\backmatter{%
  \@mainmatterfalse}

% This is for the math department, after all.
\RequirePackage{amsmath,amssymb, amsthm}
\RequirePackage{setspace}
\RequirePackage[nottoc,section]{tocbibind}

% One inch margins, one column text, USpaper
\setlength\oddsidemargin{0pt}
\setlength\textwidth{6.5in}
\setlength\topmargin{0pt}
\setlength\headheight{0pt}
\setlength\headsep{0pt}
\setlength\topskip{10pt}
% This gets the page number where it needs to be, because
% we want \texheight + 28.8pt (two baselines to the bottom of the
% page number) to give 9.25 in. I.e., page number is .75in from edge.
% this should be 8.85 (1/8 in larger) but somehow my printer puts the
% page number too low with that value.
\setlength\textheight{8.85in}
\setlength\footskip{28.8pt}

% No page headers

% Upright numbers in enumerates, even in theorems. For more
% flexibility in number formating, load the enumerate package.
\renewcommand\labelenumi{\normalfont\theenumi.}
\renewcommand\labelenumii{\normalfont(\theenumii)}
\renewcommand\labelenumiii{\normalfont\theenumiii.}
\renewcommand\labelenumiv{\normalfont\theenumiv.}

% Required: single-spaced entries, with double spacing between. Thus,
% turn off double spacing and give \itemsep value
% 2017: apparently new requirement that the references appear in the TOC

\renewenvironment{thebibliography}[1]
     {\clearpage
      \singlespacing
      \section*{\refname}%
      \addcontentsline{toc}{section}{References}
      \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
      \setlength\itemsep{.75\baselineskip plus .3333\baselineskip minus
      .1667\baselineskip}
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}

  \theoremstyle{plain}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{lemma}[theorem]{Lemma}
  \theoremstyle{definition}
\newtheorem{definition}{Definition}[section]
\newtheorem{example}{Example}[section]
  \theoremstyle{remark}
\newtheorem{remark}{Remark}


\def\bachelorinstitutionone#1{\def\@bachelorinstitutionone{#1}}
\def\bachelordegreeone#1{\def\@bachelordegreeone{#1}}
\def\bacheloryearone#1{\def\@bacheloryearone{#1}}

\def\bachelorinstitutiontwo#1{\def\@bachelorinstitutiontwo{#1}}
\def\bachelordegreetwo#1{\def\@bachelordegreetwo{#1}}
\def\bacheloryeartwo#1{\def\@bacheloryeartwo{#1}}

\def\bachelorinstitutionthree#1{\def\@bachelorinstitutionthree{#1}}
\def\bachelordegreethree#1{\def\@bachelordegreethree{#1}}
\def\bacheloryearthree#1{\def\@bacheloryearthree{#1}}

\def\masterinstitutionone#1{\def\@masterinstitutionone{#1}}
\def\masterinstitutiontwo#1{\def\@masterinstitutiontwo{#1}}
\def\masterinstitutionthree#1{\def\@masterinstitutionthree{#1}}

\def\masterdegreeone#1{\def\@masterdegreeone{#1}}
\def\masterdegreetwo#1{\def\@masterdegreetwo{#1}}
\def\masterdegreethree#1{\def\@masterdegreethree{#1}}

\def\masteryearone#1{\def\@masteryearone{#1}}
\def\masteryeartwo#1{\def\@masteryeartwo{#1}}
\def\masteryearthree#1{\def\@masteryearthree{#1}}

%Compatability
\let\bachelordegree=\bachelordegreeone
\let\bachelorinstitution=\bachelorinstitutionone
\let\bacheloryear=\bacheloryearone
\let\masterdegree=\masterdegreeone
\let\masterinstitution=\masterinstitutionone
\let\masteryear=\masteryearone

\def\dissertationdirectorone#1{\def\@dissertationdirectorone{#1}}
%\def\dissertationdirectortwo#1{\def\@dissertationdirectortwo{#1}}
\def\committeememberone#1{\def\@committeememberone{#1}}
\def\committeemembertwo#1{\def\@committeemembertwo{#1}}
%\def\committeememberthree#1{%
%   \def\@committeememberthree{#1}%
%   \def\extracommittee##1##2{##1}%
%}
%\def\committeememberfour#1{%
%   \def\@committeememberfour{#1}%
%   \def\extracommittee##1##2{##1##2}%
%}
%initialization
\let\@title\@empty
\let\@author\@empty
\let\@bachelorinstitutionone\@empty
\let\@bachelordegreeone\@empty
\let\@bacheloryearone\@empty
\let\@bachelorinstitutiontwo\@empty
\let\@bachelordegreetwo\@empty
\let\@bacheloryeartwo\@empty
\let\@bachelorinstitutionthree\@empty
\let\@bachelordegreethree\@empty
\let\@bacheloryearthree\@empty

\let\@masterinstitutionone\@empty
\let\@masteryearone\@empty
\let\@masterdegreeone\@empty
\let\@masterinstitutiontwo\@empty
\let\@masteryeartwo\@empty
\let\@masterdegreetwo\@empty
\let\@masterinstitutionthree\@empty
\let\@masteryearthree\@empty
\let\@masterdegreethree\@empty

\let\@dissertationdirector\@empty
\let\@committeememberone\@empty
\let\@committeemembertwo\@empty
\let\@committeememberthree\@empty
\let\@committeememberfour\@empty
\let\@date\@empty

\def\maketitle{%
\pagestyle{empty}
\begingroup
\clearpage
\singlespacing
\begin{centering}
\@title

\nobreak
\vspace{.55in minus .15in}
A dissertation submitted in partial fulfillment\\*
of the requirements for the degree of\\*
Master of Science in Statistics and Analytics\par

\nobreak
\vspace{.55in minus .15in}
by


\nobreak
\vspace{.55in}% minus .15in}
\@author\\*
%\vspace{0.2in}
\@bachelorinstitutionone\\*
\@bachelordegreeone, \@bacheloryearone
% more than one bachelor's degree?
\ifx\@bachelordegreetwo\@empty\else\\*
\@bachelorinstitutiontwo\\*
\@bachelordegreetwo, \@bacheloryeartwo
\fi
\ifx\@bachelordegreethree\@empty\else\\*
\@bachelorinstitutionthree\\*
\@bachelordegreethree, \@bacheloryearthree
\fi
% Up to three master's degrees
\ifx\@masterdegreeone\@empty\else\\*
\@masterinstitutionone\\*
\@masterdegreeone, \@masteryearone
\fi
\ifx\@masterdegreetwo\@empty\else\\*
\@masterinstitutiontwo\\*
\@masterdegreetwo, \@masteryeartwo
\fi
\ifx\@masterdegreethree\@empty\else\\*
\@masterinstitutionthree\\*
\@masterdegreethree, \@masteryearthree
\fi

\nobreak
\vspace{.55in minus .15in}
\@date\\*
University of Arkansas

\nobreak
\end{centering}

\nobreak
\vspace{.55in minus .15in}
\noindent This dissertation is approved for recommendation to the
Graduate Council.

\vspace{0.2in}

\nobreak
\vspace{.5in minus .1in}
\noindent \rule{3in}{.4pt}\hfil\break
\hbox to 3.24in{ \@dissertationdirectorone\hfil}\hfil\break
\hbox to 3.24in{ Thesis Director   \hfil}


\vspace{0.3in}

\nobreak
\vspace{.5in minus .1in}
\vbox{
\noindent \rule{3in}{.4pt}\hspace{.24in}\rule{3in}{.4pt}\hfil\break
\hbox to 3.24in{ \@committeememberone\hfil}%
\hbox to 3.24in{ \@committeemembertwo\hfil}\hfil\break
\hbox to 3.24in{ Committee Member   \hfil}%
\hbox to 3.24in{ Committee Member   \hfil}}%\par

\ifx\@committeememberthree\@empty\else
\nobreak
\vspace{.5in minus.1in}
\vbox{
  \extracommittee%
    {\noindent \rule{3in}{.4pt}}%
    {\hspace{.24in}\rule{3in}{.4pt}}\hfil\break
  \extracommittee%
    {\hbox to 3.24in{ \@committeememberthree\hfil}}%
    {\hbox to 3.24in{ \@committeememberfour\hfil}}\hfil\break
  \extracommittee%
    {\hbox to 3.24in{ Committee Member   \hfil}}%
    {\hbox to 3.24in{ Committee Member   \hfil}}\par
}\fi

\endgroup
\let\@title\@empty
\let\@author\@empty
%
\let\@bachelorinstitutionone\@empty
\let\@bachelordegreeone\@empty
\let\@bacheloryearone\@empty
\let\@bachelorinstitutiontwo\@empty
\let\@bachelordegreetwo\@empty
\let\@bacheloryeartwo\@empty
\let\@bachelorinstitutionthree\@empty
\let\@bachelordegreethree\@empty
\let\@bacheloryearthree\@empty
%
\let\@masterinstitutionone\@empty
\let\@masterinstitutiontwo\@empty
\let\@masterinstitutionthree\@empty
\let\@masteryearone\@empty
\let\@masteryeartwo\@empty
\let\@masteryearthree\@empty
\let\@masterdegreeone\@empty
\let\@masterdegreetwo\@empty
\let\@masterdegreethree\@empty
%
\let\@dissertationdirectorone\@empty
\let\@dissertationdirectortwo\@empty
\let\@committeememberone\@empty
\let\@committeemembertwo\@empty
\let\@committeememberthree\@empty
\let\@committeememberfour\@empty
\let\@date\@empty

\clearpage
\doublespacing
}

\renewenvironment{abstract}%
{\clearpage
 {\noindent \textbf{Abstract} \par}}
{\newpage}

\newenvironment{acknowledgements}%
{\clearpage
{\noindent \textbf{Acknowledgements} \par}}
{\newpage}

\renewcommand\tableofcontents{%
    \clearpage
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }

\raggedright
\parindent 20pt

\endinput
