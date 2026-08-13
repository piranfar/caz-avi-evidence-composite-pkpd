import zipfile, re, sys, html
def para_text(p):
    # keep tabs as separators
    out=[]
    for m in re.finditer(r'<w:(t|tab|br)(?: [^>]*)?(?:/>|>(.*?)</w:\1>)', p, re.S):
        tag=m.group(1)
        if tag=='t': out.append(html.unescape(m.group(2) or ''))
        elif tag=='tab': out.append('\t')
        elif tag=='br': out.append('\n')
    return ''.join(out)

def extract(path):
    z=zipfile.ZipFile(path)
    xml=z.read('word/document.xml').decode('utf8')
    lines=[]
    # walk body: tables become rows
    for m in re.finditer(r'<w:tbl[ >].*?</w:tbl>|<w:p[ >].*?</w:p>|<w:p/>', xml, re.S):
        blk=m.group(0)
        if blk.startswith('<w:tbl'):
            lines.append('[TABLE]')
            for r in re.finditer(r'<w:tr[ >].*?</w:tr>', blk, re.S):
                cells=[para_text(c) for c in re.findall(r'<w:tc[ >].*?</w:tc>', r.group(0), re.S)]
                lines.append(' | '.join(x.strip().replace('\n',' ') for x in cells))
            lines.append('[/TABLE]')
        else:
            lines.append(para_text(blk))
    return '\n'.join(lines)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    print(extract(sys.argv[1]))
