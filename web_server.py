import io
import os
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyramid.config import Configurator
from pyramid.renderers import render_to_response
from pyramid.response import FileResponse
from waitress import serve

from converter import common


def main_page(request):
    context = {}
    try:
        if request.POST:
            weight = float(request.POST['weight'])
            # file_name = request.POST['input-file'].filename
            input_file = request.POST['input-file']
            input_buf = io.StringIO(input_file.file.read().decode('utf-8'))
            seq_raw, seq_acids, data = common.get_data(input_buf)
            better_seq = common.get_analysis(seq_acids, weight, data)
            record = SeqRecord(better_seq, id=seq_raw.id, description=seq_raw.description)
            temp_path = os.path.join(os.getcwd(), 'temp')
            if not os.path.exists(temp_path):
                os.mkdir(temp_path)
            temp_output = tempfile.NamedTemporaryFile(dir=temp_path)
            # Release control of the file right away so that we can write to it using other libraries
            temp_output.close()
            output_path = temp_output.name
            SeqIO.write(record, output_path, 'fasta')
            response = FileResponse(output_path, request=request)
            response.headers['Content-Disposition'] = (
                f'attachment; filename={input_file.filename.replace(".fasta", "")}_optimized.fasta')
            response.headers['Content-Type'] = 'text/plain; charset'
            return response
    except:
        print('An error occurred!')
        context[ 'error'] = 'An error occurred, sorry. Check that the weight is a decimal between 0.0 and 1.0, and that the uploaded file is in the fasta format'
    return render_to_response('templates/dashboard.jinja2', context, request)


if __name__ == '__main__':
    with Configurator() as config:
        config.add_route('dashboard', '/')
        config.include('pyramid_jinja2')
        config.add_view(main_page, route_name='dashboard')
        config.add_static_view(name='static', path='templates/static')
        app = config.make_wsgi_app()
    serve(app, host='localhost', port=8888)
