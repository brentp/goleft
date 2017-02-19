package indexcov

const chartTemplate = `<!DOCTYPE html>
<html>
    <head>
{{ $name := index . "name" }}
{{ $has_sex := index . "hasSex" }}
	<title>{{ $name }}:indexcov</title>
		<script src="{{ index . "JQuery" }}"></script>
		<script src="{{ index . "ChartJS" }}"></script>
		<style type="text/css">
section {
    width: 96%;
    height: 425px;
    margin: auto;
    padding: 8px;
}

.help {
	font-family: Lucida Console;
	font-face: bold;
	font-size: +2em;
    padding: 2px;
	border: 1px solid #aaa
}


.tt {
	font-family: Lucida Console;
    border: 2px solid #aaa;
    padding: 2px;
	height: auto;
}

.one {
    width: 48%;
    height: 380px;
    float: left;
    padding: 2px;
}
.two {
    width: 48%;
    margin-left: 48%;
    height: 380px;
    padding: 2px;
}
		</style>

    </head>
    <body>
<span class="top-help">
This is the <a href="https://github.com/brentp/goleft/tree/master/indexcov">indexcov</a> summary page created with version {{ index . "version" }}.<br/>
Click the <span class="help">?</span> above each plot for help describing that type of plot.
</span>


	<section>
	<div class="one">
	<span class="tt">Inferred sex</span> <a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-sex.md" target="_blank">?</a>
	{{ if $has_sex }}
	<canvas id="canvas-sex" style="height:380px;width:380px"></canvas>
	{{ else }}
	<p> Fewer than 2 sex chromosomes found; sex plot not shown.</p>
	{{ end }}
	</div>

	<div class="two">
	<span class="tt">Bin Counts</span> <a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-bin.md" target="_blank">?</a>
	<canvas id="canvas-bin" style="height:380px;width:380px"></canvas>
	</div>

	</section>
	<hr>

<!-- links to download files; need to override height -->
<section style="height:auto">
	<div class="one" style="height:auto">
	<span class="tt">Pedigree File</span>
	<p>contains inferred sex, bins counts, and PCA values used to make the above plots</p>
	<a href="{{ $name }}-indexcov.ped">{{ $name }}-indexcov.ped</a>
	</div>

	<div class="two" style="height:auto">
	<span class="tt">Coverage BED File</span>
	<p>contains scaled coverage for every sample (each column) for each 16,384 interval in the index</p>
	<a href="{{ $name }}-indexcov.bed.gz">{{ $name }}-indexcov.bed.gz</a>
	</div>

</section><hr/>


{{ if index . "hasPCA" }}

<section style="height:auto">
	<div class="one">
	<span class="tt">PCA: 1 vs 2</span>
	<a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-pca.md" target="_blank">?</a>
	<canvas id="canvas-pca" style="height:380px;width:380px"></canvas>
	</div>

	<div class="two">
	<span class="tt">PCA: 1 vs 3</span>
	<a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-pca.md" target="_blank">?</a>
	<canvas id="canvas-pcb" style="height:380px;width:380px"></canvas>
	</div>

</section><hr/>

{{ end }}

<section style="height:auto">
{{ $notmany := index . "notmany" }}

	<div class="one">
	<span class="tt">Coverage Plots</span> <a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-depth.md#coverage" target="_blank">?</a>
	<p>Click each plot for an interactive view.</p>

	{{ $chroms := index . "chroms" }}
	{{ range $idx, $chrom := $chroms }}
		<p>
		<a href="{{ $name }}-indexcov-roc-{{ $chrom }}.html"><img src="{{ $name }}-indexcov-roc-{{ $chrom }}.png" /></a>
		</p>
	{{ end }}

<hr/>
<h5>Acknowledgements</h5>
<ul>
	<li>created with <a href="https://github.com/brentp/goleft">goleft indexcov (version {{ index . "version" }})</a> in <a href="https://golang.org">the go programming language</a></li>
	<li>bam index parsing with <a href="https://github.com/biogo/hts">biogo/hts</a></li>
	<li>interactive plots use <a href="http://www.chartjs.org/">chartjs(2)</a> via <a href="https://github.com/brentp/go-chartjs">go-chartjs</a></li>
	<li>static plots and PCA by <a href="https://github.com/gonum/plot">gonum/plot</a>and <a href="https://github.com/gonum/matrix">gonum/matrix</a></li>
</ul>

	</div>
	<div class="two">
	<span class="tt">Depth Plots</span> <a class="help" href="https://github.com/brentp/goleft/blob/master/docs/indexcov/help-depth.md" target="_blank">?</a>
	<p>
	{{ if $notmany }} Click each plot for an interactive view. {{ else }} . {{ end }}
	</p>

	{{ range $idx, $chrom := $chroms }}
		<p>
{{ if $notmany }}
		<a href="{{ $name }}-indexcov-depth-{{ $chrom }}.html"><img src="{{ $name }}-indexcov-depth-{{ $chrom }}.png" /></a>
{{ else }}
		<img src="{{ $name }}-indexcov-depth-{{ $chrom }}.png" />
{{ end }}

		</p>
	{{ end }}
	</div>


    </body>
    <script>
	Chart.defaults.line.cubicInterpolationMode = 'monotone';
	Chart.defaults.global.animation.duration = 0;

    {{ $has_sex := index . "hasSex" }}

{{ if $has_sex }}
    {{ $sex_json := index . "sex" }}
	var sex_ctx = document.getElementById("canvas-sex").getContext("2d");
	var sex_chart = new Chart(sex_ctx, {{ $sex_json }});
	var chart = sex_chart
	{{ index . "sexjs" }}
{{ end }}


    {{ $bin_json := index . "bin" }}
	var bin_ctx = document.getElementById("canvas-bin").getContext("2d");
	var bin_chart = new Chart(bin_ctx, {{ $bin_json }});
	var chart = bin_chart
	{{ index . "binjs" }}

    {{ $pca_json := index . "pca" }}
	var pca_ctx = document.getElementById("canvas-pca").getContext("2d");
	var pca_chart = new Chart(pca_ctx, {{ $pca_json }});
	var chart = pca_chart
	{{ index . "pcajs" }}

    {{ $pcb_json := index . "pcb" }}
	var pcb_ctx = document.getElementById("canvas-pcb").getContext("2d");
	var pcb_chart = new Chart(pcb_ctx, {{ $pcb_json }});
	var chart = pcb_chart
	{{ index . "pcbjs" }}

    </script>
</html>
`
