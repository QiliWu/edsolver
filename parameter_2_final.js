function crystal_system(a, b, c, alpha, beta, gamma) {
    
        if ((alpha == 90) && (beta == 90) && (gamma == 90)) {
            if ((a == b) && (b == c)) {
                return 'cubic';
            }
            if (a != b && b != c && a != c) {
                return 'orthorhombic';
            } else {
                return 'tetragonal';
            }
        }
        else if (a == b && alpha == beta && beta == 90 && gamma == 120) {
            return 'hexagonal';
        }
        else if (a != b && b != c && a != c && alpha != beta && beta != gamma && alpha != gamma) {
            return 'triclinic';
        }
        else if (a != b && b != c && a != c && alpha == gamma == 90 != beta) {
            return 'monoclinic';
        } else {
            swal({
                title:"",
                text:"晶胞参数没有对应的晶系！",
                type: "warning"
            })
            return;
        }
    
    }
    
    function d_caculate(a, b, c, alpha, beta, gamma, h, k, l) {
        var crys_sys = crystal_system(a, b, c, alpha, beta, gamma); 
        try{
            // console.log(crys_sys);
            var alpha = alpha * Math.PI / 180;
            var beta = beta * Math.PI / 180;
            var gamma = gamma * Math.PI / 180;
            // console.log(alpha,beta,gamma);
            var square_g = 0;
            if (crys_sys == 'cubic') {
                square_g = (h * h + k * k + l * l) / (a * a);
                
            }
            else if (crys_sys == 'orthorhombic') {
                square_g = (h / a) ** 2 + (k / b) ** 2 + (l / c) ** 2;
                
            }
            else if (crys_sys == 'tetragonal') {
                square_g = (h * h + k * k) / (a * a) + (l / c) ** 2;
               
            }
            else if (crys_sys == 'hexagonal') {
                square_g = 4 * (h * h + h * k + k * k) / (3 * a * a) + (l / c) ** 2;
               
            }
        
            else if (crys_sys == 'triclinic') {
                var square_v = ((a * b * c) ** 2) * (1 - Math.cos(alpha) * Math.cos(alpha) - Math.cos(beta) * Math.cos(beta) -
                    Math.cos(gamma) * Math.cos(gamma) + 2 * Math.cos(alpha) * Math.cos(beta) * Math.cos(gamma));
                var square_g = (1 / square_v) *
                    ((h * b * c * Math.sin(alpha)) ** 2 + (k * a * c * Math.sin(beta)) ** 2 + (l * b * a * Math.sin(gamma)) ** 2
                        + 2 * h * k * a * b * (c ** 2) * (Math.cos(alpha) * Math.cos(beta) - Math.cos(gamma))
                        + 2 * k * l * b * c * (a ** 2) * (Math.cos(gamma) * Math.cos(beta) - Math.cos(alpha))
                        + 2 * l * h * a * c * (b ** 2) * (Math.cos(alpha) * Math.cos(gamma) - Math.cos(beta)));
               
            }
            // console.log(square_g);
            else if (crys_sys == 'monoclinic') {
                var square_g = (h / (a * Math.sin(beta))) ** 2 + (k / b) ** 2 + (l / (c * Math.sin(beta))) ** 2 -
                    2 * h * l * Math.cos(beta) / (a * c * (Math.sin(beta)) ** 2);
                    
            }
            return (1 / Math.sqrt(square_g)).toFixed(3);
        } catch(e) {
        
            return;
        }
    }
    
    //a3
    function theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2) {
        var crys_sys = crystal_system(a, b, c, alpha, beta, gamma);
        try {
            var d_1 = d_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1);
        var d_2 = d_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2);
        // console.log(d_1,d_2);
        var alpha = alpha * Math.PI / 180;
        var beta = beta * Math.PI / 180;
        var gamma = gamma * Math.PI / 180;
        var cos_theta = 0;
        if (crys_sys == 'cubic') {
            cos_theta = (h_1 * h_2 + k_1 * k_2 + l_1 * l_2) * d_1 * d_2 / (a ** 2);
        }
    
        else if (crys_sys == 'orthorhombic') {
            cos_theta = (h_1 * h_2 / (a * a) + k_1 * k_2 / (b * b) + l_1 * l_2 / (c * c)) * d_1 * d_2;
        }
    
        else if (crys_sys == 'tetragonal') {
            cos_theta = ((h_1 * h_2 + k_1 * k_2) / (a * a) + l_1 * l_2 / (c * c)) * d_1 * d_2;
        }
    
        else if (crys_sys == 'hexagonal') {
            cos_theta = (4 / (3 * a * a)) * (h_1 * h_2 + k_1 * k_2 + 0.5 * (h_1 * k_2 + h_2 * k_1) + l_1 * l_2 / (c * c)) * d_1 * d_2;
        }
        else if (crys_sys == 'triclinic') {
            var square_v = ((a * b * c) ** 2) * (1 - Math.cos(alpha) * Math.cos(alpha) - Math.cos(beta) * Math.cos(beta)
                - Math.cos(gamma) * Math.cos(gamma) + 2 * Math.cos(alpha) * Math.cos(beta) * Math.cos(gamma));
            var A = (h_1 * h_2 * (b * c * Math.sin(alpha)) ** 2 + k_1 * k_2 * (a * c * Math.sin(beta)) ** 2 + l_1 * l_2 * (a * b * Math.sin(gamma)) ** 2
                + a * b * c * c * (Math.cos(alpha) * Math.cos(beta) - Math.cos(gamma)) * (k_1 * h_2 + h_1 * k_2)
                + a * b * b * c * (Math.cos(alpha) * Math.cos(gamma) - Math.cos(beta)) * (h_1 * l_2 + h_2 * l_1)
                + a * a * b * c * (Math.cos(beta) * Math.cos(gamma) - Math.cos(alpha)) * (k_1 * l_2 + k_2 * l_1)
            ) / square_v;
            cos_theta = A * d_1 * d_2;
        }
        else if (crys_sys == 'monoclinic') {
            cos_theta = (h_1 * h_2 / ((a * Math.sin(beta)) ** 2) + k_1 * k_2 / (b * b) + (l_1 * l_2 / ((c * Math.sin(beta)) ** 2))
                - (h_2 * l_1 + h_1 * l_2) * Math.cos(beta) / (a * c * (Math.sin(beta) ** 2)))
                * d_1 * d_2;
        }
        return (Math.acos(cos_theta) * 180 / Math.PI).toFixed(2);
        } catch(e) {
            return;
        }
        
    }
    
    //a4
    function uvw_caculate(h_1, k_1, l_1, h_2, k_2, l_2) {
        u = k_1 * l_2 - k_2 * l_1;
        v = l_1 * h_2 - l_2 * h_1;
        w = h_1 * k_2 - k_1 * h_2;
        var max_value = Math.max(Math.abs(u), Math.abs(v), Math.abs(w))
        // var arr = [u,v,w];
        // var arr1 = [];
        // arr.forEach(function(arr2){
        //     arr1.push(Math.abs(arr2));
        // })
        // console.log(arr1);
        for (let i = max_value; i > 1; i--) {
            if (u % i == 0 && v % i == 0 && w % i == 0) {
                var u = u / i, v = v / i, w = w / i;
            }
            return [parseInt(u), parseInt(v), parseInt(w)];
        }
    }
    
    function third_spot_info(h_1, k_1, l_1, h_2, k_2, l_2) {
        h_3 = h_2 - h_1;
        k_3 = k_2 - k_1;
        l_3 = l_2 - l_1;
        return [h_3, k_3, l_3];
    }
    
    // console.log(uvw_caculate(1, 1, 1, 0, 0, 2));
    // console.log(third_spot_info(1, 1, 1, 0, 0, 2));
    
    //a5
    function d_to_hkl(a, b, c, alpha, beta, gamma, d, delta_d) {
        var up = d + delta_d;
        var low = d - delta_d;
        var result = [];
        var h = 1;
        var k = 1;
        var l = 1;
        try {
            while (d_caculate(a, b, c, alpha, beta, gamma, h, 0, 0) >= low) {
                h = h + 1;
            }
            while (d_caculate(a, b, c, alpha, beta, gamma, 0, k, 0) >= low) {
                k = k + 1;
            }
            while (d_caculate(a, b, c, alpha, beta, gamma, 0, 0, l) >= low) {
                l = l + 1;
            }
            // console.log(h, k, l);
        
        
            for (let i = h - 1; i > -h; i--) {
                for (let j = k - 1; j > -k; j--) {
                    for (let p = l - 1; p > -l; p--) {
                        if ((i == j) && (j == p) && (p == 0)) {
                            continue;
                        }
                        if (d_caculate(a, b, c, alpha, beta, gamma, i, j, p) >= low &&
                            d_caculate(a, b, c, alpha, beta, gamma, i, j, p) <= up) {
                            result.push([i, j, p]);
                        }
                        // console.log(result);
                    }
                }
            }
            if (result) {
                return result;
            } else {
                swal({
                    title:"",
                    text:"找不到d值对应晶面，请尝试增加Δd!",
                    type: "warning"
                })
                return;
            }
        } catch (e) {
            return;
        }
    }
    
    //a6
    function find_hkl(a, b, c, alpha, beta, gamma, d_1, d_2, delta_d, theta, delta_theta) {
        var result = [];
        var hkl_1_list = d_to_hkl(a, b, c, alpha, beta, gamma, d_1, delta_d);
       
    
        if (hkl_1_list ) {
            var hkl_2_list = d_to_hkl(a, b, c, alpha, beta, gamma, d_2, delta_d);
            if (hkl_2_list) {
                for (s of hkl_1_list) {
                    var h_1 = s[0], k_1 = s[1], l_1 = s[2];
                    for (p of hkl_2_list) {
                        var h_2 = p[0], k_2 = p[1], l_2 = p[2];
                        // console.log(h_1,k_1,l_1,h_2,k_2,l_2);
                        if (h_1 * k_2 == h_2 * k_1 && h_1 * l_2 == h_2 * l_1 && k_1 * l_2 == k_2 * l_1) {
                            continue;
                        }
                        if ([[-h_1, -k_1, -l_1], [-h_2, -k_2, -l_2]] in result) {
                            continue;
                        } else {
                            if (theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2) <=
                                theta + delta_theta && theta_caculate(a, b, c, alpha, beta, gamma,
                                    h_1, k_1, l_1, h_2, k_2, l_2) >= theta - delta_theta) {
                                result.push([[h_1, k_1, l_1], [h_2, k_2, l_2]]);
                            }
                        }
                    }
                }
                if (result) {
                    return result;
                } else {
                    swal({
                        title:"",
                        text:"找不到匹配的h1k1l1和h2k2l2组合，请尝试增加Δd和Δθ",
                        type: "warning"
                    })
                    return;
                }
            } else {
                return;
            }
        }
    }

    // 检查输入值的类型
    function type_check(args){
        for (i in args) {
            if (isNaN(Number(args[i])) || args[i] =="") {
                swal({
                    title:"",
                    text:"请输入数字！",
                    type: "warning"
                })
                return;
            }
            if ( Number(args[i]) <= 0) {
                swal({
                    title:"",
                    text:"请输入正数！",
                    type: "warning"
                })
                return;
            }
        }
        return true;
    }   

    // 添加角度值范围的限定函数
    function angle_check(theta_1, d_1, d_2, delta_d, delta_theta) {
        if (theta_1 > 90) {
            swal({
                title:"",
                text:"请输入一个锐角的θ1",
                type: "warning"
            })
            return;
        }
        if (delta_d/d_1 > 0.3 || delta_d/d_2 > 0.3) {
            swal({
                title:"",
                text:"d值允许误差已超过30%, 结果不准确！",
                type: "warning"
            })
            return;
        }
        if (delta_theta/theta_1 > 0.3) {
            swal({
                title:"",
                text:"θ值允许误差已超过30%, 结果不准确！",
                type: "warning"
            })
            return;
        }
        return true;
    }
    // find_hkl(a, b, c, alpha, beta, gamma, d_1, d_2, delta_d, theta, delta_theta)

  //a7
function main(a, b, c, alpha, beta, gamma, d_1, d_2, delta_d, theta_1, theta_2, delta_theta) {
    if (crystal_system(a,b,c,alpha,beta,gamma)){
    var args = {a: a, b: b, c: c, alpha: alpha, beta: beta, gamma: gamma, d_1: d_1, d_2: d_2, delta_d: delta_d, theta_1: theta_1, theta_2: theta_2, delta_theta: delta_theta};
    if (type_check(args) && angle_check(theta_1, d_1, d_2, delta_d, delta_theta)){
        var hkl_ls = find_hkl(a, b, c, alpha, beta, gamma, d_1, d_2, delta_d, theta_1, delta_theta);
        var result = {};
        var arr_d1 = [];
        var arr_d2 = [];
        var arr_d3 = [];
        if (hkl_ls) {
            var n = 1;
            for (i of hkl_ls) {
                //console.log(i);)
                var h_1=i[0][0], k_1=i[0][1], l_1=i[0][2], h_2=i[1][0], k_2=i[1][1], l_2 =i[1][2]; 
               // console.log(h_1);
                var d_1_theoretical = d_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1);
                var d_2_theoretical = d_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2);
                //console.log(d_1_theoretical);
                var theta_1_theoretical = theta_caculate(a, b, c, alpha, beta, gamma, h_1, k_1, l_1, h_2, k_2, l_2);
               // console.log(theta_1_theoretical);
                var delta_theta_1 = (Math.abs(theta_1_theoretical - theta_1)).toFixed(2);
                var uvw = uvw_caculate(h_1, k_1, l_1, h_2, k_2, l_2);
                var hkl_3 = third_spot_info(h_1, k_1, l_1, h_2, k_2, l_2);
                //console.log(hkl_3);
                var d_3_theoretical= d_caculate(a, b, c, alpha, beta, gamma, hkl_3[0], hkl_3[1], hkl_3[2]);
                var theta_2_theoretical = theta_caculate(a, b, c, alpha, beta, gamma, h_2, k_2, l_2, hkl_3[0], hkl_3[1], hkl_3[2]);
                //console.log(theta_2_theoretical);
                if (Math.abs(theta_2_theoretical - theta_2) <= delta_theta) {
                   // console.log('Hello World')
                   if((arr_d1.includes(d_1_theoretical))&&(arr_d2.includes(d_2_theoretical))&&(arr_d3.includes(d_3_theoretical))){
                       continue;
                   }else{
                      result[n] = {
                        'd1': d_1_theoretical, 'd2': d_2_theoretical, 'd3': d_3_theoretical,
                        'h1k1l1': [h_1, k_1, l_1], 'h2k2l2': [h_2, k_2, l_2], 'h3k3l3': hkl_3,
                        'uvw': uvw, 'theta_1': theta_1_theoretical, 'theta_2': theta_2_theoretical,
                        'delta_theta_1': delta_theta_1
                    };
                   // console.log(result);
                   arr_d1[arr_d1.length] = d_1_theoretical;
                   arr_d2[arr_d2.length] = d_2_theoretical;
                   arr_d3[arr_d3.length] = d_3_theoretical;
                   n = n + 1;
                }
            }
            }
            if (result && result['1']){
                return result;
            } else {
                swal({
                    title:"",
                    text:"对不起，未能计算出匹配的结果",
                    type: "warning"
                })
            }
            
        } 
    }
}
    else {
        return;
    }
}

