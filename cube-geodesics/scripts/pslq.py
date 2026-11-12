import mpmath as mp
mp.mp.dps = 42
P4 = mp.mpf('0.155880815524256750509896159877085822690382')
Vcut = mp.mpf('0.00024520791519792715178055295863327505236447')
IB = mp.mpf('0.0386914290892751795440427317981596517820409303458246948129916')
s2, s3, s7 = mp.sqrt(2), mp.sqrt(3), mp.sqrt(7)
Ti2 = lambda z: mp.im(mp.polylog(2, 1j*z))
bases = {
 'elementary-A': {'1': 1, 'pi': mp.pi, 'pi*sqrt3': mp.pi*s3, 'log2': mp.log(2), 'log3': mp.log(3), 'sqrt2': s2, 'sqrt3': s3, 'sqrt7': s7},
 'elementary-B': {'1': 1, 'pi': mp.pi, 'pi*sqrt3': mp.pi*s3, 'log2': mp.log(2), 'log3': mp.log(3), 'log7': mp.log(7), 'sqrt7*atan(sqrt7)': s7*mp.atan(s7), 'sqrt7*atan(1/sqrt7)': s7*mp.atan(1/s7), 'atan(1/2)': mp.atan(mp.mpf(1)/2)},
 'dilog-A': {'1': 1, 'pi': mp.pi, 'pi*sqrt3': mp.pi*s3, 'log2': mp.log(2), 'log3': mp.log(3), 'pi^2': mp.pi**2, 'log2^2': mp.log(2)**2, 'log3^2': mp.log(3)**2, 'log2*log3': mp.log(2)*mp.log(3), 'Catalan': mp.catalan, 'Cl2(pi/3)': mp.clsin(2, mp.pi/3)},
 'dilog-B': {'1': 1, 'pi': mp.pi, 'log2': mp.log(2), 'log3': mp.log(3), 'pi*sqrt3': mp.pi*s3, 'Catalan': mp.catalan, 'pi*log2': mp.pi*mp.log(2), 'pi*log3': mp.pi*mp.log(3), 'Ti2(1/sqrt3)': Ti2(1/s3), 'Ti2(1/2)': Ti2(mp.mpf(1)/2), 'Ti2(1/3)': Ti2(mp.mpf(1)/3), 'atan(1/2)': mp.atan(mp.mpf(1)/2)},
}
for target, val in (('P4', P4), ('Vcut', Vcut), ('I_B', IB)):
    for bn, b in bases.items():
        names = list(b)
        rel = mp.pslq([val]+[mp.mpf(b[n]) for n in names], maxcoeff=2000, maxsteps=200000)
        if rel and rel[0] != 0:
            norm = mp.sqrt(sum(r*r for r in rel)); expected = mp.mpf(10)**(mp.mpf(40)/(len(names)+1))
            print(f"{target:5s} {bn:13s}: relation with norm {float(norm):.0f} (chance level for 40 digits ~{float(expected):.0f}):", dict(zip([target]+names, rel)))
        else:
            print(f"{target:5s} {bn:13s}: no relation with coefficients <= 2000")
