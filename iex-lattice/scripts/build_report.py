import json
R = json.load(open('results.json'))
def bars(vals, ses, labels, W=720, H=300, ymin=-0.45, ymax=0.25, unit='₹/kWh'):
    L, Rm, T, B = 46, 10, 14, 34
    pw, ph = W-L-Rm, H-T-B
    y = lambda v: T + (ymax-v)/(ymax-ymin)*ph
    n = len(vals); bw = pw/n
    s = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Bar chart" class="chart">']
    t = ymin
    ticks = [-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2]
    for tk in ticks:
        if abs(tk*10 - round(tk*10)) > 1e-6: continue
        s.append(f'<line x1="{L}" x2="{W-Rm}" y1="{y(tk):.1f}" y2="{y(tk):.1f}" class="{"zero" if abs(tk)<1e-9 else "grid"}"/>')
        s.append(f'<text x="{L-6}" y="{y(tk)+4:.1f}" text-anchor="end" class="tick">{tk:+.1f}</text>')
    for i,(v,se,lab) in enumerate(zip(vals,ses,labels)):
        x = L + i*bw + 2; w = bw-4
        y0, y1 = sorted([y(0), y(v)])
        s.append(f'<g class="bar"><rect x="{x:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{max(y1-y0,1):.1f}" rx="2" class="{"pos" if v>0 else "neg"}"><title>{lab}: {v:+.3f} {unit} (±{1.96*se:.3f})</title></rect>')
        cx = x+w/2
        s.append(f'<line x1="{cx:.1f}" x2="{cx:.1f}" y1="{y(v-1.96*se):.1f}" y2="{y(v+1.96*se):.1f}" class="err"/></g>')
        if i % 2 == 0: s.append(f'<text x="{cx:.1f}" y="{H-12}" text-anchor="middle" class="tick">{lab}</text>')
    s.append('</svg>'); return '\n'.join(s)
hours = [f'{h:02d}:00' for h in range(1,24)]
chart = bars(R['RTM']['byhour'], R['RTM']['se'], hours)
html = open('scripts/template.html').read().replace('{{CHART_HOUR}}', chart)
open('report.html','w').write(html)
print(len(html))
