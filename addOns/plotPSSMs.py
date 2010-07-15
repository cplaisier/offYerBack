from weblogolib import *
from weblogolib import ColorScheme
from weblogolib import ColorGroup
import numpy, corebio, cPickle, os
from pssm import pssm

def convertPSSM(pssmMatrix):
    output = ['\tPOSITION\tA\tC\tG\tT']
    for i in range(1,(len(pssmMatrix)+1)):
        output.append(str(i)+'\t'+str(pssmMatrix[i-1][0])+'\t'+str(pssmMatrix[i-1][1])+'\t'+str(pssmMatrix[i-1][2])+'\t'+str(pssmMatrix[i-1][3]))
    return output

def sortFloat(list1):
    for i in range(len(list1)):
        list1[i] = float(list1[i])
    list1.sort()
    for i in range(len(list1)):
        list1[i] = str(list1[i])
    return list1

#for file in os.listdir('randPSSMs'):
for file in ['pssms_upstream_30.pkl']:
    print file
    inFile = open('randPSSMs/'+file,'r')
    pssms = cPickle.load(inFile)
    inFile.close()
    dirSet = 'pssmPlots/'+file.rstrip('.pkl')
    if not os.path.exists(dirSet):
        os.mkdir(dirSet)
    htmlFile = open(dirSet+'/index.html','w')
    htmlFile.write('<html>\n<head><title>'+file+'</title>\n<style type=\'text/css\'>\n.hidden {\ndisplay: none;\n }\n</style\n</head>\n<body><table border=2><tr><th>Name</th><th>E-Value</th><th>PSSM</th><th>Concensus</th></tr>\n')
    pssmNames = pssms.keys()
    pssmNames = sortFloat(pssmNames)
    appendEnd = ''
    for i in pssmNames:
        if float(pssms[i].getEValue()) <= float(10):
            print i
            dist = numpy.array( pssms[i].getMatrix(), numpy.float64 ) 
            data = LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist*100)
            options = LogoOptions()
            options.Title = 'A Logo Title'
            options.color_scheme = colorscheme.nucleotide
            format = LogoFormat(data, options)
            fout = open(dirSet+'/'+str(i)+'.png', 'w')
            png_formatter(data, format, fout)
            fout.close()
            htmlFile.write('<tr><td>'+str(i)+'</td>')
            htmlFile.write('<td>')
            if float(pssms[i].getEValue()) <= float(10):
                htmlFile.write('<font color=\'#FF0000\'><b>'+str(float(pssms[i].getEValue()))+'</b></font></td>')
            else:
                htmlFile.write(str(float(pssms[i].getEValue()))+'</td>')
            htmlFile.write('<td><img src=\''+str(i)+'.png\'></td>')
            htmlFile.write('<td>'+str(pssms[i].getConsensusMotif())+'</td>')
            htmlFile.write('</tr>\n')
            appendEnd += '<div id=\''+str(i)+' style=\'display:none;\' class=\'gaggle-data ratios\'>\n'
            appendEnd += '<span class=\'gaggle-name hidden\'>'+str(i)+'</span>\n'
            appendEnd += '<span class=\'gaggle-species hidden\'>Homo sapiens</span>\n'
            appendEnd += '<span class=\'gaggle-size hidden\'>'+str(len(pssms[i].getMatrix()))+'x4</span>\n'
            appendEnd += '<div class=\'gaggle-matrix-tsv hidden\'>\n'
            appendEnd += '\n'.join(convertPSSM(pssms[i].getMatrix()))
            appendEnd += '</div>\n</div>\n'
    htmlFile.write('</table>')
    htmlFile.write(appendEnd+'</body>\n</html>')
    htmlFile.close()

#fin = open('cap.fa')
#seqs = read_seq_data(fin)
#data = LogoData.from_seqs(seqs)
#nucs = ColorScheme([ColorGroup("G", "orange"), ColorGroup("TU", "red"), ColorGroup("C", "blue"), ColorGroup("A", "green")])
#options.color_scheme = std_color_schemes["classic"]
#fout = open('caps.eps', 'w')
#eps_formatter(data, format, fout)
#fout.close()
#fout = open('caps.pdf', 'w')
#pdf_formatter(data, format, fout)
#fout.close()
#fin.close()

