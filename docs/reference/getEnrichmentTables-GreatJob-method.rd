<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Get enrichment tables from GREAT web server — getEnrichmentTables-GreatJob-method • rGREAT</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Get enrichment tables from GREAT web server — getEnrichmentTables-GreatJob-method"><meta property="og:description" content="Get enrichment tables from GREAT web server"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">
    

    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">rGREAT</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">2.1.8</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li>
  <a href="../reference/index.html">Reference</a>
</li>
<li class="dropdown">
  <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
    Articles
     
    <span class="caret"></span>
  </a>
  <ul class="dropdown-menu" role="menu"><li>
      <a href="../articles/local-GREAT.html">Analyze with local GREAT</a>
    </li>
    <li>
      <a href="../articles/online-GREAT.html">Analyze with online GREAT</a>
    </li>
    <li>
      <a href="../articles/other-docs.html">Other documents</a>
    </li>
    <li>
      <a href="../articles/other-geneset-databases.html">Work with other geneset databases</a>
    </li>
    <li>
      <a href="../articles/other-organisms.html">Work with other organisms</a>
    </li>
  </ul></li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/jokergoo/rGREAT/" class="external-link">
    <span class="fab fa-github fa-lg"></span>
     
  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->

      

      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Get enrichment tables from GREAT web server</h1>
    
    <div class="hidden name"><code>getEnrichmentTables-GreatJob-method.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Get enrichment tables from GREAT web server</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="co"># S4 method for GreatJob</span></span>
<span><span class="fu"><a href="getEnrichmentTables-dispatch.html">getEnrichmentTables</a></span><span class="op">(</span><span class="va">object</span>, ontology <span class="op">=</span> <span class="cn">NULL</span>, category <span class="op">=</span> <span class="st">"GO"</span>,</span>
<span>    request_interval <span class="op">=</span> <span class="fl">10</span>, max_tries <span class="op">=</span> <span class="fl">100</span>, download_by <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="st">"json"</span>, <span class="st">"tsv"</span><span class="op">)</span>,</span>
<span>    verbose <span class="op">=</span> <span class="cn">TRUE</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <dl><dt>object</dt>
<dd><p>A <code><a href="GreatJob-class.html">GreatJob-class</a></code> object returned by <code><a href="submitGreatJob.html">submitGreatJob</a></code>.</p></dd>

  <dt>ontology</dt>
<dd><p>Ontology names. Valid values are in <code><a href="availableOntologies-GreatJob-method.html">availableOntologies</a></code>. <code>ontology</code> is prior to  <code>category</code> argument.</p></dd>

  <dt>category</dt>
<dd><p>Pre-defined ontology categories. One category can contain more than one ontologies. Valid values are in  <code><a href="availableCategories-GreatJob-method.html">availableCategories</a></code></p></dd>

  <dt>request_interval</dt>
<dd><p>Time interval for two requests. Default is 300 seconds.</p></dd>

  <dt>max_tries</dt>
<dd><p>Maximal times for automatically reconnecting GREAT web server.</p></dd>

  <dt>download_by</dt>
<dd><p>Internally used. The complete enrichment table is provided as json data on the website, but there is no information of gene-region association. By setting <code>download_by = 'tsv'</code>, another URL from GREAT will be envoked which also contains detailed information of which genes are associated with each input region, but due to the size of the output, only top 500 terms will be returned. So if you do not really want the gene-region association column, take the default value of this argument. The columns that contain statistics are identical.</p></dd>

  <dt>verbose</dt>
<dd><p>Whether to print messages.</p></dd>


</dl></div>
    <div id="value">
    <h2>Value</h2>
    

<p>The structure of the data frames are same as the tables available on GREAT website.</p>
    </div>
    <div id="see">
    <h2>See</h2>
    <p><code><a href="availableOntologies-GreatJob-method.html">availableOntologies</a></code>, <code><a href="availableCategories-GreatJob-method.html">availableCategories</a></code></p>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>Zuguang gu &lt;z.gu@dkfz.de&gt;</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="va">job</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/readRDS.html" class="external-link">readRDS</a></span><span class="op">(</span><span class="fu"><a href="https://rdrr.io/r/base/system.file.html" class="external-link">system.file</a></span><span class="op">(</span><span class="st">"extdata"</span>, <span class="st">"GreatJob.rds"</span>, package <span class="op">=</span> <span class="st">"rGREAT"</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">tbl</span> <span class="op">=</span> <span class="fu"><a href="getEnrichmentTables-dispatch.html">getEnrichmentTables</a></span><span class="op">(</span><span class="va">job</span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The default enrichment table does not contain informatin of associated</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> genes for each input region. You can set `download_by = 'tsv'` to</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> download the complete table, but note only the top 500 regions can be</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> retreived. See the following link:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655401/Export#Export-GlobalExport</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Except the additional gene-region association column if taking 'tsv' as</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> the source of result, all other columns are the same if you choose</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> 'json' (the default) as the source. Or you can try the local GREAT</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> analysis with the function `great()`.</span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/base/names.html" class="external-link">names</a></span><span class="op">(</span><span class="va">tbl</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> [1] "GO Molecular Function" "GO Biological Process" "GO Cellular Component"</span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/utils/head.html" class="external-link">head</a></span><span class="op">(</span><span class="va">tbl</span><span class="op">[[</span><span class="fl">1</span><span class="op">]</span><span class="op">]</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>           ID</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1 GO:0070696</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2 GO:0033612</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3 GO:0070700</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4 GO:0039706</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5 GO:0043997</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6 GO:0016628</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>                                                                                    name</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1                        transmembrane receptor protein serine/threonine kinase binding</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2                                              receptor serine/threonine kinase binding</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3                                                                  BMP receptor binding</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4                                                                   co-receptor binding</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5                                  histone acetyltransferase activity (H4-K12 specific)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6 oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Binom_Genome_Fraction Binom_Expected Binom_Observed_Region_Hits</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1          3.455733e-03    3.455733000                         11</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2          3.596413e-03    3.596413000                         11</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3          2.138904e-03    2.138904000                          8</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4          2.267472e-03    2.267472000                          8</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5          3.014312e-06    0.003014312                          1</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6          2.076881e-03    2.076881000                          7</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Binom_Fold_Enrichment Binom_Region_Set_Coverage Binom_Raw_PValue</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1              3.183116                     0.011     0.0008981541</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2              3.058603                     0.011     0.0012301660</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3              3.740233                     0.008     0.0016393810</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4              3.528158                     0.008     0.0023418690</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5            331.750600                     0.001     0.0030097780</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6              3.370439                     0.007     0.0054752450</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Binom_Adjp_BH Hyper_Total_Genes Hyper_Expected Hyper_Observed_Gene_Hits</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1             1                13     0.99380020                        5</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2             1                15     1.14669300                        5</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3             1                 9     0.68801550                        4</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4             1                10     0.76446170                        4</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5             1                 1     0.07644617                        1</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6             1                24     1.83470800                        5</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Hyper_Fold_Enrichment Hyper_Gene_Set_Coverage Hyper_Term_Gene_Coverage</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1              5.031192            0.0035260930                0.3846154</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2              4.360367            0.0035260930                0.3333333</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3              5.813822            0.0028208740                0.4444444</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4              5.232440            0.0028208740                0.4000000</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5             13.081100            0.0007052186                1.0000000</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6              2.725229            0.0035260930                0.2083333</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Hyper_Raw_PValue Hyper_Adjp_BH</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 1      0.001982437     0.4919942</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2      0.004066220     0.6862153</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3      0.003134656     0.6297673</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 4      0.004910112     0.7832638</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 5      0.076446170     1.0000000</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 6      0.032459480     1.0000000</span>
<span class="r-in"><span><span class="va">job</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Submit time: 2023-04-01 09:44:07 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Note the results may only be avaiable on GREAT server for 24 hours.</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Version: 4.0.4 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Species: hg19 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Inputs: 1000 regions</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Mode: Basal plus extension </span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   Proximal: 5 kb upstream, 1 kb downstream,</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   plus Distal: up to 1000 kb</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Include curated regulatory domains</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Enrichment tables for following ontologies have been downloaded:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   GO Biological Process</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   GO Cellular Component</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>   GO Molecular Function</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">tbl</span> <span class="op">=</span> <span class="fu"><a href="getEnrichmentTables-dispatch.html">getEnrichmentTables</a></span><span class="op">(</span><span class="va">job</span>, ontology <span class="op">=</span> <span class="st">"GO Molecular Function"</span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The default enrichment table does not contain informatin of associated</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> genes for each input region. You can set `download_by = 'tsv'` to</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> download the complete table, but note only the top 500 regions can be</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> retreived. See the following link:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655401/Export#Export-GlobalExport</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Except the additional gene-region association column if taking 'tsv' as</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> the source of result, all other columns are the same if you choose</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> 'json' (the default) as the source. Or you can try the local GREAT</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> analysis with the function `great()`.</span>
<span class="r-in"><span><span class="va">tbl</span> <span class="op">=</span> <span class="fu"><a href="getEnrichmentTables-dispatch.html">getEnrichmentTables</a></span><span class="op">(</span><span class="va">job</span>, category <span class="op">=</span> <span class="st">"GO"</span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The default enrichment table does not contain informatin of associated</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> genes for each input region. You can set `download_by = 'tsv'` to</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> download the complete table, but note only the top 500 regions can be</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> retreived. See the following link:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655401/Export#Export-GlobalExport</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Except the additional gene-region association column if taking 'tsv' as</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> the source of result, all other columns are the same if you choose</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> 'json' (the default) as the source. Or you can try the local GREAT</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> analysis with the function `great()`.</span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by Zuguang Gu.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.0.7.</p>
</div>

      </footer></div>

  


  

  </body></html>

