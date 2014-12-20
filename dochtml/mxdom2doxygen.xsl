<?xml version="1.0" encoding="utf-8"?>

<!--
This is an XSL stylesheet which converts mscript XML files into XSLT.
Use the XSLT command to perform the conversion.

Ned Gulley and Matthew Simoneau, September 2003
Copyright 1984-2007 The MathWorks, Inc. 
$Revision: 1.1.6.7 $  $Date: 2008/02/29 12:45:30 $
-->

<!DOCTYPE xsl:stylesheet [ <!ENTITY nbsp "&#160;"> ]>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:escape="http://www.mathworks.com/namespace/latex/escape"
  xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd">
  <xsl:output method="text" indent="no"/>

<xsl:variable name="dTitle" select="//steptitle[@style='document']"/>
  
<xsl:template match="mscript">/*! \page <xsl:value-of select="m-file"/><xsl:text> </xsl:text><!--xsl:if test="$dTitle"><xsl:value-of select="$dTitle"/></xsl:if-->
    <!-- Determine if the there should be an introduction section. -->
    <xsl:variable name="hasIntro" select="count(cell[@style = 'overview'])"/>
    <xsl:if test = "$hasIntro">
\section Overview<xsl:text>
</xsl:text><xsl:apply-templates select="cell[1]/steptitle"/>
<xsl:apply-templates select="cell[1]/text"/>
</xsl:if>
    <xsl:variable name="body-cells" select="cell[not(@style = 'overview')]"/>
    <!-- Loop over each cell -->
    <xsl:for-each select="$body-cells">
<xsl:choose>
<xsl:when test="steptitle != 'See also'"> 
\section<xsl:text> </xsl:text><xsl:apply-templates select="steptitle"/><xsl:text>
</xsl:text><xsl:apply-templates select="text|mcode|mcodeoutput|img"/>
</xsl:when>
<xsl:otherwise>
\sa<xsl:text>
</xsl:text>
<xsl:for-each select="(text|mcode)//text()">
<xsl:call-template name="tokenize">
<xsl:with-param name="text" select="."/>
</xsl:call-template><xsl:text> </xsl:text>
</xsl:for-each>
</xsl:otherwise>
</xsl:choose>

        <!-- Contents of each cell -->
    </xsl:for-each>
*/
</xsl:template>

<xsl:template name="tokenize">
<xsl:param name="text"/>
<xsl:if test="$text!=''">
<xsl:variable name="post" select="substring-after($text,', ')"/>
<xsl:choose>
<xsl:when test="$post != ''">\ref <xsl:value-of select="substring-before($text,', ')"/>, <xsl:call-template name="tokenize">
<xsl:with-param name="text" select="$post"/>
</xsl:call-template></xsl:when>
<xsl:otherwise>\ref <xsl:value-of select="$text"/></xsl:otherwise>
</xsl:choose>
</xsl:if>
</xsl:template>


<xsl:template name="contents">
  <xsl:param name="body-cells"/>
\subsection*{Contents}

\begin{itemize}
\setlength{\itemsep}{-1ex}<xsl:for-each select="$body-cells">
      <xsl:if test="./steptitle">
   \item <xsl:apply-templates select="steptitle"/>
      </xsl:if>
    </xsl:for-each>
\end{itemize}
</xsl:template>




<!-- HTML Tags in text sections -->
<xsl:template match="p"><xsl:text>
</xsl:text>
<xsl:apply-templates/><xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="ul">\begin{itemize}
\setlength{\itemsep}{-1ex}
<xsl:apply-templates/>\end{itemize}
</xsl:template>
<xsl:template match="ol">\begin{enumerate}
\setlength{\itemsep}{-1ex}
<xsl:apply-templates/>\end{enumerate}
</xsl:template>
<xsl:template match="li">   * <xsl:apply-templates/><xsl:text>
</xsl:text></xsl:template>

<xsl:template match="pre">
  <xsl:choose>
    <xsl:when test="@class='error'">
\verbatim<xsl:text>
</xsl:text><xsl:value-of select="."/>\endverbatim
    </xsl:when>
    <xsl:otherwise>
\verbatim<xsl:text>
</xsl:text><xsl:value-of select="."/>\endverbatim
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>
<xsl:template match="b">*<xsl:apply-templates/>*</xsl:template>
<xsl:template match="tt">`<xsl:apply-templates/>`</xsl:template>
<xsl:template match="i">/<xsl:apply-templates/>/</xsl:template>
<xsl:template match="a">\verbatim<xsl:text>
</xsl:text><xsl:value-of select="."/>\endverbatim</xsl:template>

<xsl:template match="text()">
  <!-- Escape special characters in text -->
  <xsl:call-template name="replace">
    <xsl:with-param name="string" select="."/>
  </xsl:call-template>
</xsl:template>

<xsl:template match="equation">
<xsl:value-of select="."/>
</xsl:template>

<xsl:template match="latex">
    <xsl:value-of select="@text" disable-output-escaping="yes"/>
</xsl:template>
<xsl:template match="html"/>


<!-- Code input and output -->

<xsl:template match="mcode">\verbatim<xsl:text>
</xsl:text><xsl:value-of select="."/>\endverbatim
</xsl:template>


<xsl:template match="mcodeoutput">
\verbatim<xsl:value-of select="."/>\endverbatim
</xsl:template>


<!-- Figure and model snapshots -->

<xsl:template match="img">
\includegraphics [width=4in]{<xsl:value-of select="@src"/>}
</xsl:template>

<!-- Colors for syntax-highlighted input code -->

<xsl:template match="mwsh:code">\verbatim<xsl:text>
</xsl:text><xsl:apply-templates/>\endverbatim
</xsl:template>
<xsl:template match="mwsh:keywords">
  <span class="keyword"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:strings">
  <span class="string"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:comments">
  <span class="comment"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:unterminated_strings">
  <span class="untermstring"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:system_commands">
  <span class="syscmd"><xsl:value-of select="."/></span>
</xsl:template>


<!-- Used to escape special characters in the LaTeX output. -->

<escape:replacements>
  <!-- special TeX characters -->
  <replace><from>$</from><to>\$</to></replace>
  <replace><from>@</from><to>\@</to></replace>
  <replace><from>&amp;</from><to>\&amp;</to></replace>
  <replace><from>%</from><to>\%</to></replace>
  <replace><from>#</from><to>\#</to></replace>
  <replace><from>_</from><to>\_</to></replace>
  <replace><from>{</from><to>\{</to></replace>
  <replace><from>}</from><to>\}</to></replace>
  <!-- mainly in code -->
  <replace><from>\</from><to>\\</to></replace>
  <!-- mainly in math -->
  <replace><from>&lt;</from><to>\&lt;</to></replace>
  <replace><from>&gt;</from><to>\&gt;</to></replace>
</escape:replacements>

<xsl:variable name="replacements" select="document('')/xsl:stylesheet/escape:replacements/replace"/>

<xsl:template name="replace">
  <xsl:param name="string"/>
  <xsl:param name="next" select="1"/>

  <xsl:variable name="count" select="count($replacements)"/>
  <xsl:variable name="first" select="$replacements[$next]"/>
  <xsl:choose>
    <xsl:when test="$next > $count">
      <xsl:value-of select="$string"/>
    </xsl:when>
    <xsl:when test="contains($string, $first/from)">      
      <xsl:call-template name="replace">
        <xsl:with-param name="string"
                        select="substring-before($string, $first/from)"/>
        <xsl:with-param name="next" select="$next+1" />
      </xsl:call-template>
      <xsl:copy-of select="$first/to" />
      <xsl:call-template name="replace">
        <xsl:with-param name="string"
                        select="substring-after($string, $first/from)"/>
        <xsl:with-param name="next" select="$next"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="replace">
        <xsl:with-param name="string" select="$string"/>
        <xsl:with-param name="next" select="$next+1"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
