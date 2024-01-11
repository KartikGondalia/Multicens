document.addEventListener('DOMContentLoaded', () => {
  var form = document.getElementById('input-form');
  var uploadButton = document.getElementById('submit');
  var outputDiv = document.getElementById('output');

  var runSampleButton = document.getElementById('sample-run');

  runSampleButton.onclick = function(event){
    var source = document.getElementById('source_tissue');
    var target = document.getElementById('target_tissue');
    var hormone = document.getElementById('hormone');
    var input1 = document.getElementById('input1');
    var input2 = document.getElementById('input2');

    source.selectedIndex = 10;
    target.selectedIndex = 8;
    hormone.selectedIndex = 6;

    // input1.value = 'insulin_source_genes_full.csv';
    // input2.value = 'insulin_target_genes_full.csv';

    outputDiv.innerHTML = '<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/>  Please wait! It may take 3-5 minutes ...';

    var xhr = new XMLHttpRequest();
      xhr.open('GET', '/run_sample', true);    

      // Set up a handler for when the task for the request is complete.
      xhr.onload = function () {
        if (xhr.status === 200) {
          outputDiv.innerHTML = xhr.responseText;
          $('#ranking').DataTable({'columnDefs': [
            { className: 'dt-center', targets: '_all' },
              ],
              'order': [2, 'asc'],
          });
        } else {
          outputDiv.innerHTML = 'An error occurred during the upload. Try again.';
        }
      };

      // Send the data.
      xhr.send();
  }

  form.onsubmit = function(event) {
      event.preventDefault();

      outputDiv.innerHTML = '<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/> Please wait! It may take 3-5 minutes ...';

      var file1 = document.getElementById('input1').files[0];
      var file2 = document.getElementById('input2').files[0];

      var formData = new FormData(form);

      if (!(file1.type.match('.csv') || file2.type.match('.csv'))) {
          outputDiv.innerHTML = 'CSV file is required!.';
          return;
      }

      if (file1.size >= 20000000 || file2.size >= 20000000) {
          outputDiv.innerHTML = 'You cannot upload this file because its size exceeds the maximum limit of 2 MB.';
          return;
      }

      // formData.append('input1', file1, file1.name);
      // formData.append('input2', file2, file2.name);

      var xhr = new XMLHttpRequest();
      xhr.open('POST', '/get_ranking', true);
  

      // Set up a handler for when the task for the request is complete.
      xhr.onload = function () {
        if (xhr.status === 200) {
          outputDiv.innerHTML = xhr.responseText;
          $('#ranking').DataTable({'columnDefs': [
                  { className: 'dt-center', targets: '_all' },
              ],
              'order': [2, 'asc'],
          });
        } else {
          outputDiv.innerHTML = 'An error occurred during the upload. Try again.';
        }
      };

      // Send the data.
      xhr.send(formData);
  }
});