package indexcov

const chartTemplate = `<!DOCTYPE html>
<html>
    <head>
		<script src="{{ index . "JQuery" }}"></script>
		<script src="{{ index . "ChartJS" }}"></script>
		<style type="text/css">

section {
    width: 96%;
    height: 425px;
    margin: auto;
    padding: 8px;
}


.tt {
	font-family: Lucida Console;
    border: 1px solid #aaa;
    padding: 1px;
}

.one {
    width: 48%;
    height: 380px;
    float: left;
}
.two {
    width: 48%;
    margin-left: 48%;
    height: 380px;
}
		</style>

    </head>
    <body>
	<section>
	<div class="one">
	<span class="tt">Inferred Sex</span>
	<canvas id="canvas-sex" style="height:380px;width:380px"></canvas>
	</div>

	<div class="two">
	<span class="tt">Bin Counts</span>
	<canvas id="canvas-bin" style="height:380px;width:380px"></canvas>
	</div>

	</section>

	<hr>

	<section>
	<div class="one">
	<span class="tt">PCA: 1 vs 2</span>
	<canvas id="canvas-pca" style="height:380px;width:380px"></canvas>
	</div>

	<div class="two">
	<span class="tt">PCA: 1 vs 3</span>
	<canvas id="canvas-pcb" style="height:380px;width:380px"></canvas>
	</div>

	</section>
    </body>
    <script>
	Chart.defaults.line.cubicInterpolationMode = 'monotone';
	Chart.defaults.global.animation.duration = 0;

    {{ $sex_json := index . "sex" }}
	var sex_ctx = document.getElementById("canvas-sex").getContext("2d");
	var sex_chart = new Chart(sex_ctx, {{ $sex_json }});
	var chart = sex_chart
	{{ index . "sexjs" }}

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
</html>`
